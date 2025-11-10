# 📚 Lane-Emden Documentation Index

**Welcome!** This is your guide to the new Lane-Emden scalable architecture.

## 🚀 Quick Start (Choose Your Path)

### Path 1: "Just show me what to do!" ⚡
1. Run: `./setup_lane_emden.sh`
2. Run: `make lane_emden_list`
3. Run: `make lane_emden PRESET=polytrope_n1.5_3d`
4. Done! ✨

### Path 2: "I want to understand first" 📖
1. Read: `THIS_IS_WHAT_YOU_ASKED_FOR.md` (5 min)
2. Read: `VISUAL_SUMMARY.md` (10 min)
3. Run: `./setup_lane_emden.sh`
4. Test: `make lane_emden PRESET=polytrope_n1.5_3d`

### Path 3: "I need all the details" 🔬
1. Read: `VISUAL_SUMMARY.md`
2. Read: `MIGRATION_GUIDE.md`
3. Read: `lane_emden/docs/LANE_EMDEN_ARCHITECTURE.md`
4. Run: `./setup_lane_emden.sh`

## 📄 Document Guide

### 🎯 Essential (Read These First)

| Document | Purpose | Reading Time | When to Read |
|----------|---------|--------------|--------------|
| **THIS_IS_WHAT_YOU_ASKED_FOR.md** | Executive summary | 5 min | Before you start |
| **VISUAL_SUMMARY.md** | Before/after comparison | 10 min | To understand the changes |
| **WORKFLOW_COMPARISON.md** | Old vs new workflows | 10 min | To see practical differences |

### 🔧 Implementation

| Document | Purpose | Reading Time | When to Read |
|----------|---------|--------------|--------------|
| **MIGRATION_GUIDE.md** | Step-by-step migration | 15 min | During migration |
| **lane_emden/README.md** | Quick reference | 5 min | Daily usage |

### 📚 Reference

| Document | Purpose | Reading Time | When to Read |
|----------|---------|--------------|--------------|
| **lane_emden/docs/LANE_EMDEN_ARCHITECTURE.md** | Complete technical design | 30 min | For deep understanding |
| **README_INDEX.md** | This file! | 5 min | To navigate docs |

## 🛠️ Scripts & Tools

### Setup Scripts

| Script | Purpose | Usage |
|--------|---------|-------|
| **setup_lane_emden.sh** | ⭐ One-command setup | `./setup_lane_emden.sh` |
| **migrate_to_new_architecture.sh** | Detailed migration | `./migrate_to_new_architecture.sh` |

### Management Tools

| Tool | Purpose | Usage |
|------|---------|-------|
| **lane_emden_config_manager.py** | Manage presets | `python3 lane_emden_config_manager.py list` |
| **Makefile.lane_emden** | Build targets | `make lane_emden PRESET=polytrope_n1.5_3d` |

## 📖 Reading Order by Role

### 👤 End User (Just want to run simulations)

1. `THIS_IS_WHAT_YOU_ASKED_FOR.md` - What changed
2. `VISUAL_SUMMARY.md` - How to use it
3. Run `./setup_lane_emden.sh`
4. `lane_emden/README.md` - Daily reference

### 🔬 Researcher (Want to create custom configurations)

1. `THIS_IS_WHAT_YOU_ASKED_FOR.md` - Overview
2. `VISUAL_SUMMARY.md` - Capabilities
3. `MIGRATION_GUIDE.md` - Setup process
4. `lane_emden/docs/LANE_EMDEN_ARCHITECTURE.md` - Configuration schema
5. Practice: Create custom preset

### 💻 Developer (Want to extend the system)

1. `THIS_IS_WHAT_YOU_ASKED_FOR.md` - Goals
2. `VISUAL_SUMMARY.md` - Current state
3. `lane_emden/docs/LANE_EMDEN_ARCHITECTURE.md` - Technical design
4. `MIGRATION_GUIDE.md` - Implementation details
5. Study: `lane_emden_config_manager.py` source code

## 🎓 Learning Path

### Beginner: First Time Using

```
1. Read: THIS_IS_WHAT_YOU_ASKED_FOR.md
   ↓
2. Run: ./setup_lane_emden.sh
   ↓
3. Try: make lane_emden_list
   ↓
4. Run: make lane_emden PRESET=polytrope_n1.5_3d
   ↓
5. Analyze results
```

### Intermediate: Want Custom Configurations

```
1. Review: VISUAL_SUMMARY.md
   ↓
2. Learn: WORKFLOW_COMPARISON.md
   ↓
3. Study: lane_emden/README.md
   ↓
4. Create: python3 lane_emden_config_manager.py create ...
   ↓
5. Validate: make lane_emden_validate PRESET=custom
   ↓
6. Run: make lane_emden PRESET=custom
```

### Advanced: Extending the System

```
1. Read: lane_emden/docs/LANE_EMDEN_ARCHITECTURE.md
   ↓
2. Study: lane_emden_config_manager.py source
   ↓
3. Understand: Configuration schema
   ↓
4. Design: New features
   ↓
5. Implement: Extensions
   ↓
6. Test: All presets
   ↓
7. Document: Your changes
```

## 🗂️ Document Categories

### 📋 Overview Documents
- `THIS_IS_WHAT_YOU_ASKED_FOR.md` - Executive summary
- `VISUAL_SUMMARY.md` - Visual overview
- `README_INDEX.md` - This file

### 🔄 Migration Documents
- `MIGRATION_GUIDE.md` - How to migrate
- `WORKFLOW_COMPARISON.md` - Old vs new

### 📘 Technical Documents
- `lane_emden/docs/LANE_EMDEN_ARCHITECTURE.md` - Full design
- `lane_emden/config/presets/*.json` - Preset examples

### 📗 User Guides
- `lane_emden/README.md` - Quick reference
- Inline help: `make lane_emden_help`

## 🔍 Find What You Need

### "How do I...?"

| Question | Answer Location |
|----------|----------------|
| ...run a simulation? | `VISUAL_SUMMARY.md` Quick Start |
| ...create a new preset? | `lane_emden/README.md` Creating Custom Configurations |
| ...understand the physics? | `lane_emden/docs/LANE_EMDEN_ARCHITECTURE.md` Preset Selection Guide |
| ...migrate my old setup? | `MIGRATION_GUIDE.md` Step-by-step |
| ...validate a config? | `lane_emden/README.md` Examples |
| ...change relaxation steps? | `WORKFLOW_COMPARISON.md` Command Flow |

### "What is...?"

| Term | Explanation Location |
|------|---------------------|
| Preset | `VISUAL_SUMMARY.md` New Solution |
| Polytropic index | `lane_emden/docs/LANE_EMDEN_ARCHITECTURE.md` Preset Selection Guide |
| Adiabatic index | `VISUAL_SUMMARY.md` Polytrope Physics Cheat Sheet |
| Config manager | `MIGRATION_GUIDE.md` What Changed |

### "Where is...?"

| Item | Location |
|------|----------|
| Presets | `lane_emden/config/presets/` |
| Data files | `lane_emden/data/numerical_solutions/` |
| Scripts | `lane_emden/scripts/` |
| Results | `lane_emden/results/` |
| Documentation | `lane_emden/docs/` |
| Old results | `lane_emden/results/legacy/` |

## 📊 Document Map

```
Root Documentation
├── README_INDEX.md (this file)            ← You are here
├── THIS_IS_WHAT_YOU_ASKED_FOR.md         ← Start here
├── VISUAL_SUMMARY.md                      ← Overview
├── WORKFLOW_COMPARISON.md                 ← Old vs new
├── MIGRATION_GUIDE.md                     ← How to migrate
│
└── lane_emden/
    ├── README.md                          ← Quick reference
    │
    ├── config/
    │   ├── presets/                       ← Example configs
    │   │   ├── polytrope_n1.5_3d.json
    │   │   └── polytrope_n1.5_2d.json
    │   └── templates/                     ← Base templates
    │       ├── base_2d.json
    │       └── base_3d.json
    │
    └── docs/
        ├── LANE_EMDEN_ARCHITECTURE.md     ← Complete design
        ├── LANE_EMDEN_NOTES.md            ← Legacy notes
        └── [other docs]
```

## 🎯 By Task

### Task: Run Your First Simulation

**Documents to read**:
1. `THIS_IS_WHAT_YOU_ASKED_FOR.md` (Quick Start section)
2. `VISUAL_SUMMARY.md` (Quick Start section)

**Commands to run**:
```bash
./setup_lane_emden.sh
make lane_emden_list
make lane_emden PRESET=polytrope_n1.5_3d
```

**Time needed**: 10 minutes

### Task: Create Custom Configuration

**Documents to read**:
1. `VISUAL_SUMMARY.md` (Preset Selection Guide)
2. `lane_emden/README.md` (Creating Custom Configurations)
3. `lane_emden/docs/LANE_EMDEN_ARCHITECTURE.md` (Configuration Schema)

**Commands to run**:
```bash
python3 lane_emden_config_manager.py create \
  --name custom_n3_3d --polytrope 3.0 --dimension 3
make lane_emden_validate PRESET=custom_n3_3d
make lane_emden PRESET=custom_n3_3d
```

**Time needed**: 30 minutes

### Task: Understand the Architecture

**Documents to read**:
1. `VISUAL_SUMMARY.md` (complete)
2. `WORKFLOW_COMPARISON.md` (complete)
3. `lane_emden/docs/LANE_EMDEN_ARCHITECTURE.md` (complete)

**Time needed**: 1 hour

### Task: Migrate Existing Setup

**Documents to read**:
1. `MIGRATION_GUIDE.md` (complete)
2. `VISUAL_SUMMARY.md` (Backward Compatibility section)

**Commands to run**:
```bash
./setup_lane_emden.sh
make lane_emden_list
make lane_emden PRESET=polytrope_n1.5_3d  # Test
```

**Time needed**: 20 minutes

## 💡 Tips

### For Quick Learners
- Start with `VISUAL_SUMMARY.md`
- Run `./setup_lane_emden.sh`
- Experiment with commands
- Refer back to docs as needed

### For Thorough Learners
- Read all overview documents first
- Understand the architecture before setup
- Study examples in detail
- Plan your custom configurations

### For Troubleshooting
- Check `MIGRATION_GUIDE.md` Troubleshooting section
- Review `lane_emden/README.md` Examples
- Validate with `make lane_emden_validate`
- Check inline help: `make lane_emden_help`

## 🔗 Quick Links

### Most Important

- 🎯 **Start Here**: `THIS_IS_WHAT_YOU_ASKED_FOR.md`
- 📖 **Learn**: `VISUAL_SUMMARY.md`
- 🚀 **Setup**: `./setup_lane_emden.sh`
- 📚 **Reference**: `lane_emden/README.md`

### For Help

- ❓ **How to use**: `WORKFLOW_COMPARISON.md`
- 🔧 **How to migrate**: `MIGRATION_GUIDE.md`
- 📐 **Technical details**: `lane_emden/docs/LANE_EMDEN_ARCHITECTURE.md`
- 🆘 **Troubleshooting**: `MIGRATION_GUIDE.md` (Troubleshooting section)

## 📞 Getting Help

1. **Check the docs** (this index!)
2. **Run inline help**: `make lane_emden_help`
3. **Validate your config**: `make lane_emden_validate PRESET=yourpreset`
4. **Review examples** in `lane_emden/config/presets/`

## ✨ What's Next?

After reading this index:

1. **New Users**: Read `THIS_IS_WHAT_YOU_ASKED_FOR.md`
2. **Migrating**: Read `MIGRATION_GUIDE.md`
3. **Learning**: Read `VISUAL_SUMMARY.md`
4. **Developing**: Read `lane_emden/docs/LANE_EMDEN_ARCHITECTURE.md`

---

**Happy simulating!** 🚀
