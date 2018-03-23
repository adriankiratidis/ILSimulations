
all:
	make -C src/
	make -C bin/

clean:
	make -C src/ clean
	make -C bin/ clean
