# eQTLEnrichment — Claude skills branch

This orphan branch holds only Claude Code skills for this repo, versioned
independently of the code branches. It is checked out once, in its own worktree
(`Benchmark/eQTLEnrichment-skills/`), and never switched or merged into a code
branch.

- `skills/eqtl-benchmark-runner/` — run, extend and interpret the eQTL benchmark.

Claude loads these through the EngreitzLabAgents plugin: the local checkout's
`skills/<name>` is a git-ignored symlink to `skills/<name>` here, so every session
(any worktree, any branch, any directory) uses this one copy. Edit and commit here.
