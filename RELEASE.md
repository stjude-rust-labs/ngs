# Release

  * [ ] Update version in `Cargo.toml`.
  * [ ] Update `CHANGELOG.md` with version and publication date.
  * [ ] Run tests: `cargo test --all-features`.
  * [ ] Run linting: `cargo clippy --all-features`.
  * [ ] Run fmt: `cargo fmt --check`.
  * [ ] Stage changes: `git add Cargo.lock Cargo.toml CHANGELOG.md`.
  * [ ] Create git commit: `git commit -m "release: bumps version to v0.1.0"`.
  * [ ] Create git tag: `git tag v0.1.0`.
  * [ ] Push release: `git push && git push --tags`.
  * [ ] Publish the new crate: `cargo publish --all-features`.
  * [ ] Go to the Releases page in Github, create a Release for this tag, and
    copy the notes from the `CHANGELOG.md` file.
  * [ ] Run docker release workflow from release tag branch.
