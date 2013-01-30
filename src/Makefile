OBJS = Main.o General.o SimData.o BEAST_XML.o User_Interface.o Analysis_Settings.o
CC = g++
DEBUG = -g
CFLAGS = -Wall -c -m64 -O3 -funroll-loops $(DEBUG)
LFLAGS = -Wall -m64 $(DEBUG)

BEASTifier: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o BEASTifier

Main.o: Main.cpp
	$(CC) $(CFLAGS) Main.cpp

General.o: General.cpp General.h
	$(CC) $(CFLAGS) General.cpp

User_Interface.o: User_Interface.cpp User_Interface.h
	$(CC) $(CFLAGS) User_Interface.cpp

Analysis_Settings.o: Analysis_Settings.cpp Analysis_Settings.h
	$(CC) $(CFLAGS) Analysis_Settings.cpp

SimData.o: SimData.cpp SimData.h
	$(CC) $(CFLAGS) SimData.cpp

BEAST_XML.o: BEAST_XML.cpp BEAST_XML.h
	$(CC) $(CFLAGS) BEAST_XML.cpp

clean:
	rm -rf *.o BEASTifier