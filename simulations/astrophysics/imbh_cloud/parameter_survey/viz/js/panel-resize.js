// ============================================================
// IMBH-Cloud Visualization - Panel Resize Functionality
// Allows dragging panel corners to resize
// ============================================================

// Panels that can be resized
const RESIZABLE_PANELS = [
    'info-panel',
    'shock-panel',
    'toggle-panel',
    'sim-select',
    'legend',
    'pv-panel',
    'selection-panel',
    'colorbar-container'
];

// Minimum sizes for panels
const MIN_PANEL_SIZE = {
    'info-panel': { width: 150, height: 100 },
    'shock-panel': { width: 130, height: 60 },
    'toggle-panel': { width: 140, height: 100 },
    'sim-select': { width: 160, height: 100 },
    'legend': { width: 120, height: 80 },
    'pv-panel': { width: 200, height: 200 },
    'selection-panel': { width: 160, height: 150 },
    'colorbar-container': { width: 50, height: 100 }
};

// State for resize operation
let resizeState = {
    isResizing: false,
    panel: null,
    startX: 0,
    startY: 0,
    startWidth: 0,
    startHeight: 0,
    handle: null
};

// Initialize resize functionality
function initPanelResize() {
    RESIZABLE_PANELS.forEach(panelId => {
        const panel = document.getElementById(panelId);
        if (panel) {
            addResizeHandle(panel);
        }
    });

    // Global mouse events for resize
    document.addEventListener('mousemove', onResizeMove);
    document.addEventListener('mouseup', onResizeEnd);

    console.log('Panel resize initialized');
}

// Add resize handle to a panel
function addResizeHandle(panel) {
    // Create resize handle element
    const handle = document.createElement('div');
    handle.className = 'resize-handle';
    handle.innerHTML = '⋮⋮'; // Visual indicator
    handle.title = 'Drag to resize, double-click to reset';

    // Store panel reference
    handle.dataset.panelId = panel.id;

    // Mouse down starts resize
    handle.addEventListener('mousedown', onResizeStart);

    // Double-click resets size
    handle.addEventListener('dblclick', (e) => {
        e.preventDefault();
        e.stopPropagation();
        resetPanelSize(panel.id);
    });

    // Append to panel
    panel.appendChild(handle);

    // Make sure panel has position relative for handle positioning
    const computedStyle = window.getComputedStyle(panel);
    if (computedStyle.position === 'static') {
        panel.style.position = 'relative';
    }
}

// Start resize operation
function onResizeStart(e) {
    e.preventDefault();
    e.stopPropagation();

    const panelId = e.target.dataset.panelId;
    const panel = document.getElementById(panelId);

    if (!panel) return;

    resizeState.isResizing = true;
    resizeState.panel = panel;
    resizeState.handle = e.target;
    resizeState.startX = e.clientX;
    resizeState.startY = e.clientY;
    resizeState.startWidth = panel.offsetWidth;
    resizeState.startHeight = panel.offsetHeight;

    // Add resizing class for visual feedback
    panel.classList.add('resizing');
    document.body.style.cursor = 'nwse-resize';

    // Disable text selection during resize
    document.body.style.userSelect = 'none';
}

// Handle resize movement
function onResizeMove(e) {
    if (!resizeState.isResizing || !resizeState.panel) return;

    const panel = resizeState.panel;
    const panelId = panel.id;
    const minSize = MIN_PANEL_SIZE[panelId] || { width: 100, height: 80 };

    // Calculate new dimensions
    const deltaX = e.clientX - resizeState.startX;
    const deltaY = e.clientY - resizeState.startY;

    let newWidth = resizeState.startWidth + deltaX;
    let newHeight = resizeState.startHeight + deltaY;

    // Apply minimum constraints
    newWidth = Math.max(minSize.width, newWidth);
    newHeight = Math.max(minSize.height, newHeight);

    // Apply maximum constraints (don't exceed viewport)
    const rect = panel.getBoundingClientRect();
    const maxWidth = window.innerWidth - rect.left - 20;
    const maxHeight = window.innerHeight - rect.top - 20;
    newWidth = Math.min(maxWidth, newWidth);
    newHeight = Math.min(maxHeight, newHeight);

    // Apply new dimensions
    panel.style.width = newWidth + 'px';
    panel.style.height = newHeight + 'px';
    panel.style.maxWidth = 'none'; // Override max-width constraints

    // Special handling for PV panel - resize canvas too
    if (panelId === 'pv-panel') {
        resizePVCanvas(newWidth, newHeight);
    }

    // Special handling for colorbar
    if (panelId === 'colorbar-container') {
        resizeColorbar(newWidth, newHeight);
    }
}

// End resize operation
function onResizeEnd(e) {
    if (!resizeState.isResizing) return;

    if (resizeState.panel) {
        resizeState.panel.classList.remove('resizing');
    }

    resizeState.isResizing = false;
    resizeState.panel = null;
    resizeState.handle = null;

    document.body.style.cursor = '';
    document.body.style.userSelect = '';
}

// Resize PV canvas when panel is resized
function resizePVCanvas(panelWidth, panelHeight) {
    const canvas = document.getElementById('pv-canvas');
    if (!canvas) return;

    // Calculate available space for canvas (accounting for padding and other elements)
    const availableWidth = panelWidth - 20;  // padding
    const availableHeight = panelHeight - 120; // space for controls below

    // Minimum canvas size
    const canvasWidth = Math.max(180, availableWidth);
    const canvasHeight = Math.max(120, availableHeight);

    canvas.width = canvasWidth;
    canvas.height = canvasHeight;
    canvas.style.width = canvasWidth + 'px';
    canvas.style.height = canvasHeight + 'px';

    // Redraw PV diagram with new size
    if (STATE.snapshots && STATE.snapshots.length > 0) {
        updatePVDiagram(STATE.snapshots[STATE.currentFrame]);
    }
}

// Resize colorbar when panel is resized
function resizeColorbar(panelWidth, panelHeight) {
    const colorbar = document.getElementById('colorbar');
    const labels = document.getElementById('colorbar-labels');

    if (colorbar) {
        const newHeight = Math.max(80, panelHeight - 30);
        colorbar.style.height = newHeight + 'px';
    }
    if (labels) {
        const newHeight = Math.max(80, panelHeight - 30);
        labels.style.height = newHeight + 'px';
    }
}

// Reset panel to default size
function resetPanelSize(panelId) {
    const panel = document.getElementById(panelId);
    if (!panel) return;

    panel.style.width = '';
    panel.style.height = '';
    panel.style.maxWidth = '';

    // Reset PV canvas if needed
    if (panelId === 'pv-panel') {
        const canvas = document.getElementById('pv-canvas');
        if (canvas) {
            canvas.width = 260;
            canvas.height = 180;
            canvas.style.width = '260px';
            canvas.style.height = '180px';
            if (STATE.snapshots && STATE.snapshots.length > 0) {
                updatePVDiagram(STATE.snapshots[STATE.currentFrame]);
            }
        }
    }
}

// Reset all panels to default sizes
function resetAllPanelSizes() {
    RESIZABLE_PANELS.forEach(panelId => {
        resetPanelSize(panelId);
    });
    console.log('All panel sizes reset to defaults');
}

// Initialize when DOM is ready
document.addEventListener('DOMContentLoaded', initPanelResize);
