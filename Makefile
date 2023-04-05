CC = gcc
RM = rm -f
CFLAGS = -Wall -Wextra -pedantic -std=c11 -g
LDFLAGS = -g 
LDLIBS = -lpari -lm

SRCS = real_quadratic_orders.c fp_representations.c
OBJS = $(subst .c,.o,$(SRCS))

testing: testing.o $(OBJS)
	$(CC) $(LDFLAGS) -o testing testing.o $(OBJS) $(LDLIBS)

testing.o: testing.c testing.h
real_quadratic_orders.o: real_quadratic_orders.c real_quadratic_orders.h
fp_representations.o: fp_representations.c fp_representations.h

clean:
	$(RM) $(OBJS) testing.o