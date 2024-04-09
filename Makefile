#
#
#
CFLAG = /usr/bin/gcc -Wall -O3 -flto -g -I./include
BIN_NAME = postmajor
SRC = ./src
INSTALL_DIR = /usr/local/bin
UTILITY = $(SRC)/iirfilter.o $(SRC)/picker_wu.o $(SRC)/sac.o $(SRC)/seisdata_load.o

#
all: postmajor
#
postmajor: $(SRC)/postmajor.o $(UTILITY)
	$(CFLAG) -o $@ $(SRC)/postmajor.o $(UTILITY) -lm

# Compile rule for Object
%.o:%.c
	$(CFLAG) -c $< -o $@

#
#
install:
	@echo Installing $(BIN_NAME) to $(INSTALL_DIR)...
	@cp ./$(BIN_NAME) $(INSTALL_DIR)
	@echo Finish installing of $(BIN_NAME).

# Clean-up rules
clean:
	(cd $(SRC); rm -f *.o *.obj *% *~; cd -)

clean_bin:
	rm -f $(BIN_NAME)
