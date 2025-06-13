
build:
	cargo build --release

install:
	cp target/release/gsmm2-metric /usr/bin/

bai:
	cargo build --release
	cp target/release/gsmm2-metric /usr/bin/