# Changelog — f2py-meson-flang-windows

Provenance for this skill. Each entry records what changed, why, and which
knowledge candidates motivated it. Format: skill-forge
`templates/changelog.md`.

```yaml
- date: 2026-09-13
  action: create
  skill_name: f2py-meson-flang-windows
  scope: project
  path: .claude/skills/f2py-meson-flang-windows/
  source_candidates:
    - KC-20260913-001
    - KC-20260913-002
  derived_from: []
  changes:
    - "create skill from end-to-end verified bhmiepy build — curator suggested
      global scope, but the evidence gate failed the spot-check (flang-Windows
      toolchain half verified in a single project, no official documentation),
      so it landed project-level per R1; promote to global via the agent_skill
      repo + skill-deploy once a second project uses it or the user asks"
    - "hard rules: .bat-not-one-liner (%LIB% parse-time expansion, KC-20260913-002),
      implementation sources must be extension sources, root-level custom_target,
      ASCII-only Fortran"
    - "troubleshooting table maps all six failure modes hit during bhmiepy init
      to their fixes"
    - "templates: dev-install.cmd and root meson.build pattern"
  authorization: ""
  validation:
    frontmatter: passed
    structure: passed
    scope: passed
    provenance: passed
```

```yaml
- date: 2026-09-13
  action: update
  skill_name: f2py-meson-flang-windows
  scope: project
  path: .claude/skills/f2py-meson-flang-windows/
  source_candidates:
    - KC-20260913-001
  derived_from: []
  changes:
    - "new Hard Rule R5 (F77-clean f2py file + external bridge to F90
      modules) — f2py parser cannot handle modules/use/interfaces/kinds in
      externally-linked code; bridge also avoids the module-proc vs
      same-named F77 external symbol collision (both bind to name_)"
    - "troubleshooting rows: flang stacks large F77 locals -> exit
      0xC00000FD stack overflow with no traceback, fix with SAVE; meson-
      python editable rebuild-on-import runs plain ninja without the
      toolchain env (WinError 2 / LNK1104), rule: full dev-install rerun
       before invoking Python"
    - "description gains the 0xC00000FD trigger phrase; meson template
      comments document the bridge pattern"
  authorization: ""
  validation:
    frontmatter: passed
    structure: passed
    scope: passed
    provenance: passed
```

```yaml
- date: 2026-09-13
  action: update
  skill_name: f2py-meson-flang-windows
  scope: project
  path: .claude/skills/f2py-meson-flang-windows/
  source_candidates:
    - KC-20260913-001
  derived_from: []
  changes:
    - "new Hard Rule R6: f2py reorders generated signatures (dimension-
      owning args become optional and move after the array; positional
      calls silently misbind) — call entry points with keyword args after
      the leading arrays, check ext.__doc__ when in doubt"
    - "editable-rebuild troubleshooting row amended: fires after Fortran
      AND meson.build edits (regeneration step exits 1 without the
      toolchain env)"
  authorization: ""
  validation:
    frontmatter: passed
    structure: passed
    scope: passed
    provenance: passed
```

```yaml
- date: 2026-09-13
  action: update
  skill_name: f2py-meson-flang-windows
  scope: project
  path: .claude/skills/f2py-meson-flang-windows/
  source_candidates:
    - KC-20260913-001
  derived_from: []
  changes:
    - "new Hard Rule R7 from the code-review round: guard EVERY exported
      entry into process-fatal (stop-based) Fortran, or return ierr;
      port guard formulas comparing in float64 before int casts (int64
      wraparound at ~9.2e18 bypasses checks); preserve comparison
      polarity for NaN semantics ('if not (w > 0)' not 'if w <= 0')"
  authorization: ""
  validation:
    frontmatter: passed
    structure: passed
    scope: passed
    provenance: passed
```

```yaml
- date: 2026-09-14
  action: update
  skill_name: f2py-meson-flang-windows
  scope: project
  path: .claude/skills/f2py-meson-flang-windows/
  source_candidates:
    - KC-20260913-001
  derived_from: []
  changes:
    - "troubleshooting rows from the Phase 2 upstream-A/B round: (a) flang 23
      aborts 'runtime derived type info descriptor was not generated' on
      derived types with allocatable components as intent(out) dummies ->
      per-target fortran_args ['-mllvm', '--ignore-missing-type-desc'] (kebab
      spelling; the diagnostic's camelCase name is rejected), validate
      outputs afterwards; (b) second 0xC00000FD flavor: MB-scale automatic
      arrays vs Windows' 1 MB stack reserve -> '-Wl,/STACK:67108864' link
      arg guarded by host_machine.system(), no source edits to pristine
      upstream code; (c) meson-python editable artifacts live in
      build/<python-tag>/ — scan build/ one level deep for non-installed
      targets; regeneration also fires after 'meson configure -D...'"
    - "description gains the RTTI-abort trigger phrase and covers running
      (not only calling) legacy Fortran"
  authorization: ""
  validation:
    frontmatter: passed
    structure: passed
    scope: passed
    provenance: passed
```
