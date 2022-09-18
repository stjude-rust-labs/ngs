# Release

  * [ ] Update version in `Cargo.toml`.
  * [ ] Update `CHANGELOG.md` with version and publication date.
  * [ ] Run tests: `cargo test`.
  * [ ] Run linting: `cargo clippy`.
  * [ ] Run fmt: `cargo fmt`.
  * [ ] Stage changes: `git add Cargo.lock Cargo.toml CHANGELOG.md`.
  * [ ] Create git commit: `git commit -m "release: bumps version to v0.1.0"`.
  * [ ] Create git tag: `git tag v0.1.0`.
  * [ ] Push release: `git push && git push --tags`.
  * [ ] Publish the new crate: `cargo publish`.