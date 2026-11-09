FROM rust:1.98.1-alpine

WORKDIR /app
COPY . .
RUN cargo build --release

FROM alpine
COPY --from=0 /app/target/release/linear-programming /usr/local/bin/linear-programming

ENTRYPOINT ["linear-programming"]