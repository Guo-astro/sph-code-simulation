#!/bin/bash
# Deploy IMBH-Cloud Visualization to GitHub Pages
# Usage: ./deploy_github_pages.sh [repo-name] [username]

set -e

REPO_NAME="${1:-imbh-cloud-viz}"
USERNAME="${2:-}"

echo "=== IMBH-Cloud Visualization Deployment ==="
echo ""

# Check if gh CLI is available
if ! command -v gh &> /dev/null; then
    echo "GitHub CLI (gh) not found. Please install it:"
    echo "  brew install gh"
    echo "  gh auth login"
    exit 1
fi

# Check if authenticated
if ! gh auth status &> /dev/null; then
    echo "Please authenticate with GitHub first:"
    echo "  gh auth login"
    exit 1
fi

# Get username if not provided
if [ -z "$USERNAME" ]; then
    USERNAME=$(gh api user -q '.login')
fi

echo "Repository: $REPO_NAME"
echo "Username: $USERNAME"
echo ""

# Check if in viz directory
if [ ! -f "imbh_cloud_viz.html" ]; then
    echo "Error: Run this script from the viz/ directory"
    exit 1
fi

# Check data exists
if [ ! -d "data" ]; then
    echo "Error: data/ directory not found. Run convert_to_binary.py first."
    exit 1
fi

# Calculate total size
TOTAL_SIZE=$(du -sh . | cut -f1)
echo "Total size to deploy: $TOTAL_SIZE"
echo ""

# Create or update repo
if gh repo view "$USERNAME/$REPO_NAME" &> /dev/null; then
    echo "Repository exists, updating..."
else
    echo "Creating new repository..."
    gh repo create "$REPO_NAME" --public --description "IMBH-Cloud Tidal Disruption Visualization" || true
fi

# Initialize git if needed
if [ ! -d ".git" ]; then
    git init
    git branch -M main
fi

# Add remote
REMOTE_URL="https://github.com/$USERNAME/$REPO_NAME.git"
if git remote | grep -q origin; then
    git remote set-url origin "$REMOTE_URL"
else
    git remote add origin "$REMOTE_URL"
fi

# Add and commit
git add .
git commit -m "Update visualization" || echo "No changes to commit"

# Push
echo ""
echo "Pushing to GitHub..."
git push -u origin main --force

# Enable GitHub Pages
echo ""
echo "Enabling GitHub Pages..."
gh api -X PUT "repos/$USERNAME/$REPO_NAME/pages" \
    -f build_type=legacy \
    -f source='{"branch":"main","path":"/"}' 2>/dev/null || \
gh api -X POST "repos/$USERNAME/$REPO_NAME/pages" \
    -f build_type=legacy \
    -f source='{"branch":"main","path":"/"}' 2>/dev/null || \
echo "Pages may already be enabled"

echo ""
echo "=== Deployment Complete ==="
echo ""
echo "Your visualization will be available at:"
echo "  https://$USERNAME.github.io/$REPO_NAME/imbh_cloud_viz.html"
echo ""
echo "Note: It may take 1-2 minutes for the site to become available."
