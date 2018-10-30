# General flags
CC = gcc
CFLAGS = -Wall -std=c99 -O2 -g 
LDFLAGS = -lm

BIN = build/rareman

OBJS = build/rareman.o build/matpbm.o
HEADERS = src/matpbm.h
SOURCES = src/rareman.c src/matpbm.c 

all:		$(BIN)

$(BIN):		$(OBJS)
		$(CC) $(CFLAGS) $(OBJS) -o $(BIN) $(LDFLAGS)

build/%.o:	src/%.c $(HEADERS)
		mkdir -p build
		$(CC) $(CFLAGS) -c $< -o $@

clean:
		rm -rf build

