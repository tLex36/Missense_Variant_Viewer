const geneInput = document.getElementById('gene-input');
const transcriptSelect = document.getElementById('transcript-select');
const generateBtn = document.getElementById('generate-btn');
const lookupStatus = document.getElementById('lookup-status');
const generateStatus = document.getElementById('generate-status');
const resultsSection = document.getElementById('results');

const windowSizeInput = document.getElementById('window-size');
const windowSizeVal = document.getElementById('window-size-val');
const normalizeMissenseInput = document.getElementById('normalize-missense');
const flipGnomadInput = document.getElementById('flip-gnomad');
const intronSizeInput = document.getElementById('intron-size');

let lastGeneLookedUp = '';

// --- Theme toggle (light/dark) ---------------------------------------------------
// The initial theme is already applied by a blocking script in index.html's <head>,
// before first paint, so there's no flash of the wrong theme - this just wires up
// the button and keeps localStorage in sync with any change.

const themeToggle = document.getElementById('theme-toggle');

function applyTheme(theme) {
    document.documentElement.setAttribute('data-theme', theme);
    localStorage.setItem('theme', theme);
    if (themeToggle) {
        themeToggle.setAttribute('aria-label', theme === 'light' ? 'Switch to dark theme' : 'Switch to light theme');
    }
}

if (themeToggle) {
    applyTheme(document.documentElement.getAttribute('data-theme') || 'dark');
    themeToggle.addEventListener('click', () => {
        const current = document.documentElement.getAttribute('data-theme') === 'light' ? 'light' : 'dark';
        applyTheme(current === 'light' ? 'dark' : 'light');
    });
}

// --- Gene -> transcript lookup -------------------------------------------------

async function lookupTranscripts() {
    const gene = geneInput.value.trim().toUpperCase();
    if (!gene || gene === lastGeneLookedUp) return;

    transcriptSelect.disabled = true;
    transcriptSelect.innerHTML = '<option value="">Looking up transcripts…</option>';
    generateBtn.disabled = true;
    lookupStatus.textContent = '';
    lookupStatus.classList.remove('error');

    try {
        const resp = await fetch(`/api/transcripts?gene=${encodeURIComponent(gene)}`);
        const data = await resp.json();

        if (!resp.ok) {
            transcriptSelect.innerHTML = '<option value="">No transcripts found</option>';
            lookupStatus.textContent = data.error || 'Gene not found.';
            lookupStatus.classList.add('error');
            return;
        }

        lastGeneLookedUp = gene;
        transcriptSelect.innerHTML = '';
        data.transcripts.forEach((t, i) => {
            const opt = document.createElement('option');
            opt.value = t.id;
            opt.textContent = t.label;
            transcriptSelect.appendChild(opt);
        });
        transcriptSelect.disabled = false;
        generateBtn.disabled = false;
    } catch (err) {
        lookupStatus.textContent = 'Lookup failed: ' + err;
        lookupStatus.classList.add('error');
    }
}

geneInput.addEventListener('change', lookupTranscripts);
geneInput.addEventListener('keydown', (e) => {
    if (e.key === 'Enter') { e.preventDefault(); lookupTranscripts(); }
});

// --- Tabs -----------------------------------------------------------------------

document.querySelectorAll('.tab-btn').forEach(btn => {
    btn.addEventListener('click', () => {
        document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
        document.querySelectorAll('.tab-panel').forEach(p => p.classList.remove('active'));
        btn.classList.add('active');
        document.getElementById('tab-' + btn.dataset.tab).classList.add('active');
    });
});

// --- Generate ---------------------------------------------------------------------

windowSizeInput.addEventListener('input', () => {
    windowSizeVal.textContent = windowSizeInput.value;
});

async function generate() {
    const gene = geneInput.value.trim().toUpperCase();
    const transcriptId = transcriptSelect.value;
    if (!gene || !transcriptId) return;

    generateBtn.disabled = true;
    generateStatus.classList.remove('error');
    generateStatus.textContent = '⏱ Generating chart… this can take 1–5 minutes for a gene that has never been requested before.';
    resultsSection.classList.remove('hidden');

    const payload = {
        gene,
        transcript_id: transcriptId,
        window_size: parseInt(windowSizeInput.value, 10),
        normalize_missense: normalizeMissenseInput.checked,
        flip_gnomad: flipGnomadInput.checked,
        intron_size: parseInt(intronSizeInput.value, 10),
    };

    try {
        const resp = await fetch('/api/generate', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(payload),
        });
        const data = await resp.json();

        if (!resp.ok) {
            generateStatus.textContent = data.error || 'Generation failed.';
            generateStatus.classList.add('error');
            return;
        }

        renderResults(data);
        generateStatus.textContent = `Done — ${data.gene} / ${data.transcript_id} (${data.length} aa)`;
    } catch (err) {
        generateStatus.textContent = 'Generation failed: ' + err;
        generateStatus.classList.add('error');
    } finally {
        generateBtn.disabled = false;
    }
}

generateBtn.addEventListener('click', generate);

// Re-generate automatically when any chart control changes, once a first result exists
[windowSizeInput, normalizeMissenseInput, flipGnomadInput].forEach(el => {
    el.addEventListener('change', () => {
        if (!resultsSection.classList.contains('hidden')) generate();
    });
});
intronSizeInput.addEventListener('change', () => {
    if (!resultsSection.classList.contains('hidden')) generate();
});

function renderResults(data) {
    // Plotly chart HTML fragments (Plotly.js is already loaded globally via CDN in index.html)
    document.getElementById('missense-chart').innerHTML = data.missense_chart_html;
    executeScripts(document.getElementById('missense-chart'));

    document.getElementById('trunc-chart').innerHTML = data.trunc_chart_html;
    executeScripts(document.getElementById('trunc-chart'));

    const missenseTableEl = document.getElementById('missense-table');
    missenseTableEl.innerHTML = data.missense_table_html;
    enhanceTable(missenseTableEl);

    const truncTableEl = document.getElementById('trunc-table');
    truncTableEl.innerHTML = data.trunc_table_html;
    enhanceTable(truncTableEl);

    // The structure viewer is a full standalone HTML document - load it via srcdoc
    // so its own <script> tags execute in a sandboxed context, same as the original
    // Streamlit components.html() iframe behavior.
    document.getElementById('structure-frame').srcdoc = data.structure_html;
}

// --- Table sorting (click a column header to sort by it; click again to reverse) ------

function enhanceTable(container) {
    const table = container.querySelector('table');
    if (!table) return;
    const thead = table.querySelector('thead');
    const tbody = table.querySelector('tbody');
    if (!thead || !tbody) return;

    // Use the last header row in case pandas rendered a multi-level header
    const headerRow = thead.rows[thead.rows.length - 1];
    Array.from(headerRow.cells).forEach((th, colIndex) => {
        th.dataset.sortDir = '';
        th.addEventListener('click', () => sortTableByColumn(table, colIndex, th));
    });
}

function sortTableByColumn(table, colIndex, th) {
    const tbody = table.querySelector('tbody');
    const rows = Array.from(tbody.querySelectorAll('tr'));
    const newDir = th.dataset.sortDir === 'asc' ? 'desc' : 'asc';

    table.querySelectorAll('thead th').forEach(h => {
        h.dataset.sortDir = '';
        h.classList.remove('sorted-asc', 'sorted-desc');
    });
    th.dataset.sortDir = newDir;
    th.classList.add(newDir === 'asc' ? 'sorted-asc' : 'sorted-desc');

    const cellValue = (row) => {
        const cell = row.children[colIndex];
        return cell ? cell.textContent.trim() : '';
    };

    rows.sort((rowA, rowB) => {
        const a = cellValue(rowA);
        const b = cellValue(rowB);
        const numA = parseFloat(a);
        const numB = parseFloat(b);
        const bothNumeric = a !== '' && b !== '' && !isNaN(numA) && !isNaN(numB);
        const cmp = bothNumeric
            ? numA - numB
            : a.localeCompare(b, undefined, { numeric: true, sensitivity: 'base' });
        return newDir === 'asc' ? cmp : -cmp;
    });

    rows.forEach(row => tbody.appendChild(row));
}

// innerHTML doesn't execute <script> tags - Plotly's fig.to_html() output includes one,
// so we re-insert it manually to get it to run.
function executeScripts(container) {
    container.querySelectorAll('script').forEach(oldScript => {
        const newScript = document.createElement('script');
        Array.from(oldScript.attributes).forEach(attr => newScript.setAttribute(attr.name, attr.value));
        newScript.textContent = oldScript.textContent;
        oldScript.parentNode.replaceChild(newScript, oldScript);
    });
}