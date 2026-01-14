// ============================================================
// IMBH-Cloud Visualization - Panel Auto-Hide (Focus Mode)
// Hides UI panels for distraction-free viewing
// ============================================================

// Configuration
const AUTOHIDE_CONFIG = {
    hideDelay: 3000,        // ms before hiding when mouse leaves window
    idleDelay: 8000,        // ms of inactivity before hiding (0 to disable)
    fadeOutDuration: 400,   // ms for fade out animation
    fadeInDuration: 200,    // ms for fade in animation
    showHintDelay: 500,     // ms before showing the hint after hiding
};

// State
let panelHideTimeout = null;
let idleTimeout = null;
let panelsHidden = false;
let mouseInWindow = true;

// All panel IDs to control
const PANEL_IDS = [
    'info-panel',
    'shock-panel',
    'toggle-panel',
    'sim-select',
    'controls',
    'colorbar-container',
    'legend',
    'view-buttons',
    'pv-panel',
    'selection-panel'
];

// Initialize auto-hide functionality
function initPanelAutoHide() {
    // Create hint element
    createHiddenHint();

    // Add CSS for transitions
    addTransitionStyles();

    // Mouse leave/enter window events
    document.addEventListener('mouseleave', onMouseLeaveWindow);
    document.addEventListener('mouseenter', onMouseEnterWindow);

    // Mouse move for idle detection
    document.addEventListener('mousemove', onMouseActivity);
    document.addEventListener('mousedown', onMouseActivity);
    document.addEventListener('wheel', onMouseActivity);

    // Keyboard shortcut
    document.addEventListener('keydown', onKeyDown);

    // Start idle timer
    resetIdleTimer();

    console.log('Panel auto-hide initialized (Press H to toggle, or move mouse outside window)');
}

// Create the hint element shown when panels are hidden
function createHiddenHint() {
    const hint = document.createElement('div');
    hint.id = 'panels-hidden-hint';
    hint.innerHTML = `
        <span style="opacity: 0.7;">Panels hidden</span>
        <span style="margin: 0 8px;">•</span>
        <span>Move mouse or press <kbd>H</kbd> to show</span>
    `;
    document.body.appendChild(hint);
}

// Add CSS transition styles
function addTransitionStyles() {
    const style = document.createElement('style');
    style.textContent = `
        /* Panel transition styles */
        #info-panel, #shock-panel, #toggle-panel, #sim-select,
        #controls, #colorbar-container, #legend, #view-buttons,
        #pv-panel, #selection-panel {
            transition: opacity ${AUTOHIDE_CONFIG.fadeOutDuration}ms ease-out,
                        transform ${AUTOHIDE_CONFIG.fadeOutDuration}ms ease-out;
        }

        .panels-hidden #info-panel,
        .panels-hidden #shock-panel,
        .panels-hidden #toggle-panel,
        .panels-hidden #sim-select,
        .panels-hidden #controls,
        .panels-hidden #colorbar-container,
        .panels-hidden #legend,
        .panels-hidden #view-buttons,
        .panels-hidden #pv-panel,
        .panels-hidden #selection-panel {
            opacity: 0 !important;
            pointer-events: none !important;
            transform: scale(0.98);
        }

        /* Hint styling */
        #panels-hidden-hint {
            position: fixed;
            bottom: 20px;
            left: 50%;
            transform: translateX(-50%) translateY(20px);
            background: rgba(20, 20, 40, 0.9);
            color: #aaa;
            padding: 10px 20px;
            border-radius: 20px;
            border: 1px solid #4a4a6a;
            font-size: 13px;
            opacity: 0;
            pointer-events: none;
            transition: opacity 0.3s, transform 0.3s;
            z-index: 1000;
        }

        #panels-hidden-hint kbd {
            background: #3a3a5a;
            padding: 2px 6px;
            border-radius: 3px;
            border: 1px solid #5a5a7a;
            font-family: monospace;
            font-size: 12px;
        }

        .panels-hidden #panels-hidden-hint {
            opacity: 1;
            transform: translateX(-50%) translateY(0);
        }

        /* Loading overlay should never hide */
        #loading {
            transition: none !important;
        }
        .panels-hidden #loading {
            opacity: 1 !important;
            pointer-events: auto !important;
            transform: none !important;
        }
    `;
    document.head.appendChild(style);
}

// Hide all panels
function hidePanels(instant = false) {
    if (panelsHidden) return;

    panelsHidden = true;

    if (instant) {
        // Skip animation
        document.body.classList.add('panels-hidden-instant');
    }

    document.body.classList.add('panels-hidden');
    console.log('Panels hidden (focus mode)');
}

// Show all panels
function showPanels() {
    if (!panelsHidden) return;

    panelsHidden = false;
    document.body.classList.remove('panels-hidden');
    document.body.classList.remove('panels-hidden-instant');

    // Reset timers
    resetIdleTimer();

    console.log('Panels shown');
}

// Toggle panels visibility
function togglePanels() {
    if (panelsHidden) {
        showPanels();
    } else {
        hidePanels();
    }
}

// Mouse left the browser window
function onMouseLeaveWindow() {
    mouseInWindow = false;

    // Clear any existing timeout
    if (panelHideTimeout) {
        clearTimeout(panelHideTimeout);
    }

    // Start hide timer
    panelHideTimeout = setTimeout(() => {
        if (!mouseInWindow) {
            hidePanels();
        }
    }, AUTOHIDE_CONFIG.hideDelay);
}

// Mouse entered the browser window
function onMouseEnterWindow() {
    mouseInWindow = true;

    // Cancel hide timeout
    if (panelHideTimeout) {
        clearTimeout(panelHideTimeout);
        panelHideTimeout = null;
    }

    // Show panels if hidden
    if (panelsHidden) {
        showPanels();
    }

    // Reset idle timer
    resetIdleTimer();
}

// Any mouse activity (move, click, scroll)
function onMouseActivity() {
    // Show panels if hidden and mouse is in window
    if (panelsHidden && mouseInWindow) {
        showPanels();
    }

    // Reset idle timer
    resetIdleTimer();
}

// Reset the idle timer
function resetIdleTimer() {
    if (AUTOHIDE_CONFIG.idleDelay <= 0) return;

    if (idleTimeout) {
        clearTimeout(idleTimeout);
    }

    idleTimeout = setTimeout(() => {
        if (mouseInWindow && !panelsHidden) {
            hidePanels();
        }
    }, AUTOHIDE_CONFIG.idleDelay);
}

// Keyboard handler
function onKeyDown(e) {
    // 'H' key to toggle panels (but not when typing in input)
    if (e.key.toLowerCase() === 'h' &&
        !e.ctrlKey && !e.altKey && !e.metaKey &&
        e.target.tagName !== 'INPUT' &&
        e.target.tagName !== 'TEXTAREA' &&
        e.target.tagName !== 'SELECT') {

        e.preventDefault();
        togglePanels();
    }

    // Escape to show panels if hidden
    if (e.key === 'Escape' && panelsHidden) {
        showPanels();
    }
}

// Check if panels are currently hidden
function arePanelsHidden() {
    return panelsHidden;
}

// Temporarily disable auto-hide (e.g., during drag operations)
function disableAutoHide() {
    if (panelHideTimeout) {
        clearTimeout(panelHideTimeout);
        panelHideTimeout = null;
    }
    if (idleTimeout) {
        clearTimeout(idleTimeout);
        idleTimeout = null;
    }
}

// Re-enable auto-hide
function enableAutoHide() {
    resetIdleTimer();
}

// Initialize when DOM is ready
document.addEventListener('DOMContentLoaded', initPanelAutoHide);
