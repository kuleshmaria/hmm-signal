CC     = gcc -m64
CFLAGS = -Wall -Wextra -Wpedantic

OBJ = main.o utils.o hmm.o baum_welch.o

main: $(OBJ)
	$(CC) $(CFLAGS) -o main $(OBJ)

clean:
	rm -f main $(OBJ)
