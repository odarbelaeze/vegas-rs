FROM rust:1.90-bookworm AS builder

WORKDIR /usr/src/vegas
COPY . .
RUN cargo install --path .

FROM debian:bookworm
COPY --from=builder /usr/local/cargo/bin/vegas /usr/local/bin/vegas
CMD [ "vegas" ]
