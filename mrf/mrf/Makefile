SRC =  mrf.cpp ICM.cpp GCoptimization.cpp graph.cpp maxflow.cpp \
       MaxProdBP.cpp LinkedBlockList.cpp regions-maxprod.cpp \
       example.cpp

BIN = example

CC = g++
CPPFLAGS = -O2 -W -Wall
OBJ = $(SRC:.cpp=.o)

all: $(BIN)

$(BIN): $(OBJ)
	$(CC) -o $(BIN) $(OBJ)

clean: 
	rm -f $(BIN) $(OBJ)
