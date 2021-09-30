#
#
#
BIN_NAME = postmajor
CFLAGS = $(GLOBALFLAGS) -O2 -g -I./include -flto

UTILITY = iirfilter.o picker_wu.o

SRC = postmajor_sac.c iirfilter.c picker_wu.c

postmajor_sac: postmajor_sac.o $(UTILITY)
	$(CC) $(CFLAGS) -o $(BIN_NAME) postmajor_sac.o $(UTILITY) -lm

# Compile rule for Object
.c.o:
	$(CC) $(CFLAGS) -c $<

# Clean-up rules
clean:
	rm -f a.out core *.o *.obj *% *~

clean_bin:
	rm -f $(BIN_NAME)
