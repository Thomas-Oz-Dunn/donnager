FROM rust:1.67

WORKDIR /usr/src/donnager
COPY . .

RUN cargo install --path .
