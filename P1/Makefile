CC=g++
CFLAGS= -Wall -Wextra -std=gnu99 -g
SFILES = server.c connect.c jobs.c send.c
CFILES = client.c connect.c jobs.c fetch.c


all: client.c server.c protocol.h connect.c job.h jobs.c fetch.c send.c
	$(CC) $(CFILES) -o client $(CFLAGS)
