#
#
#
CFLAG = /usr/bin/gcc -Wall -O3 -flto -g -I./include
BIN_NAME = postmajor chan_mod
SRC = ./src
UTILITY = $(SRC)/iirfilter.o $(SRC)/picker_wu.o $(SRC)/sac_proc.o

all: postmajor chan_mod

postmajor: $(SRC)/postmajor_sac.o $(UTILITY)
	$(CFLAG) -o $@ $(SRC)/postmajor_sac.o $(UTILITY) -lm

chan_mod: $(SRC)/chan_mod.o $(SRC)/sac_proc.o
	$(CFLAG) -o $@ $(SRC)/chan_mod.o $(SRC)/sac_proc.o

concat_sac: $(SRC)/concat_sac.o $(SRC)/sac_proc.o
	$(CFLAG) -o $@ $(SRC)/concat_sac.o $(SRC)/sac_proc.o

preproc_sac: $(SRC)/preproc_sac.o $(SRC)/sac_proc.o
	$(CFLAG) -o $@ $(SRC)/preproc_sac.o $(SRC)/sac_proc.o

integral_sac: $(SRC)/integral_sac.o $(SRC)/sac_proc.o
	$(CFLAG) -o $@ $(SRC)/integral_sac.o $(SRC)/sac_proc.o

# Compile rule for Object
%.o:%.c
	$(CFLAG) -c $< -o $@

# Clean-up rules
clean:
	(cd $(SRC); rm -f *.o *.obj *% *~; cd -)

clean_bin:
	rm -f $(BIN_NAME)
