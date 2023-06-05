CC = gcc
RM = rm -f
CFLAGS = -Wall -Wextra -pedantic -std=c11 -g -fPIC
LDFLAGS = -g
LDLIBS = -lpari

SRCS = real_quadratic_orders.c fp_representations.c compact_representations.c twin_smooths.c utility.c 
OBJS = $(subst .c,.o,$(SRCS))

testing: testing.o $(OBJS)
	$(CC) $(LDFLAGS) -o testing testing.o $(OBJS) $(LDLIBS)

shared: $(OBJS)
	$(CC) $(LDFLAGS) -shared -o infrastructure.so $(OBJS) $(LDLIBS)

testing.o: testing.c testing.h
utility.o: utility.c utility.h
real_quadratic_orders.o: real_quadratic_orders.c real_quadratic_orders.h
fp_representations.o: fp_representations.c fp_representations.h
compact_representations.o: compact_representations.c compact_representations.h
twin_smooths.o: twin_smooths.c twin_smooths.h

clean:
	$(RM) $(OBJS) testing.o testing infrastructure.so
