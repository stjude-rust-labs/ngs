# This Dockerfile builds on debian:buster-slim to include the ngs tool
# from the current working directory.
#
# This produces a docker imange that contains a worker ngs binary. Images
# released on the Github container registry are build using this file.

FROM rust:1.64-slim-buster AS builder

RUN apt-get update && apt-get install -y git

WORKDIR /usr/src/ngs

COPY Cargo.toml .
COPY Cargo.lock .
COPY .cargo ./.cargo
COPY src ./src

RUN cargo install --path .

FROM debian:buster-slim
LABEL maintainer="Keivn Benton <krbenton.opensource@icloud.com>"

COPY --from=builder /usr/local/cargo/bin/ngs /usr/local/bin/ngs

WORKDIR /data

ENTRYPOINT ["ngs"]
