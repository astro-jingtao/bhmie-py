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
