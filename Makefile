#
#
#
CFLAG = /usr/bin/gcc -Wall -O3 -flto -g -I./include
BIN_NAME = postmajor
SRC = ./src/
UTILITY = iirfilter.o picker_wu.o sac_proc.o

postmajor_sac: postmajor_sac.o $(UTILITY)
	$(CFLAG) -o $(BIN_NAME) postmajor_sac.o $(UTILITY) -lm

# Compile rule for Object
%.o:$(SRC)%.c
	$(CFLAG) -c $< -o $@

# Clean-up rules
clean:
	rm -f a.out core *.o *.obj *% *~

clean_bin:
	rm -f $(BIN_NAME)
