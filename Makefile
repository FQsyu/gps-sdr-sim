# Makefile for Linux etc.

.PHONY: all clean time
all: gps-sdr-sim

SHELL=/bin/bash
CC=gcc
NVCC?=nvcc
CFLAGS=-O3 -Wall -D_FILE_OFFSET_BITS=64
ifdef USER_MOTION_SIZE
CFLAGS+=-DUSER_MOTION_SIZE=$(USER_MOTION_SIZE)
endif
LDFLAGS=-lm

# Check if CUDA is enabled and available
ENABLE_CUDA ?= $(shell which $(NVCC) >/dev/null 2>&1 && echo 1 || echo 0)

ifeq ($(ENABLE_CUDA), 1)
CFLAGS+=-DENABLE_CUDA
NVCCFLAGS=-O3 -arch=native
OBJS=gpssim.o l5_code.o cnav.o gpssim_cuda.o
LD=$(NVCC)
else
OBJS=gpssim.o l5_code.o cnav.o
LD=$(CC)
endif

gps-sdr-sim: $(OBJS)
	${LD} $^ ${LDFLAGS} -o $@

test_l5_code: test_l5_code.o l5_code.o
	${CC} $^ ${LDFLAGS} -o $@

test_l5_code.o: test_l5_code.c l5_code.h
	${CC} -c ${CFLAGS} $< -o $@

test_cnav: test_cnav.o cnav.o l5_code.o
	${CC} $^ ${LDFLAGS} -o $@

test_cnav.o: test_cnav.c cnav.h l5_code.h
	${CC} -c ${CFLAGS} $< -o $@

l5_code.o: l5_code.c l5_code.h
	${CC} -c ${CFLAGS} $< -o $@

cnav.o: cnav.c cnav.h gpssim.h
	${CC} -c ${CFLAGS} $< -o $@

gpssim.o: gpssim.c .user-motion-size gpssim.h l5_code.h cnav.h
	${CC} -c ${CFLAGS} $< -o $@

gpssim_cuda.o: gpssim_cuda.cu gpssim_cuda.cuh gpssim.h
	${NVCC} -c ${NVCCFLAGS} $< -o $@

.user-motion-size: .FORCE
	@if [ -f .user-motion-size ]; then \
		if [ "`cat .user-motion-size`" != "$(USER_MOTION_SIZE)" ]; then \
			echo "Updating .user-motion-size"; \
			echo "$(USER_MOTION_SIZE)" >| .user-motion-size; \
		fi; \
	else \
		echo "$(USER_MOTION_SIZE)" > .user-motion-size; \
	fi;

clean:
	rm -f gpssim.o l5_code.o cnav.o gpssim_cuda.o test_l5_code.o test_l5_code test_cnav.o test_cnav gps-sdr-sim *.bin .user-motion-size



time: gps-sdr-sim
	time ./gps-sdr-sim -e brdc3540.14n -u circle.csv -b 1
	time ./gps-sdr-sim -e brdc3540.14n -u circle.csv -b 8
	time ./gps-sdr-sim -e brdc3540.14n -u circle.csv -b 16

.FORCE:

YEAR?=$(shell date +"%Y")
Y=$(patsubst 20%,%,$(YEAR))
%.$(Y)n:
	wget -q ftp://cddis.gsfc.nasa.gov/gnss/data/daily/$(YEAR)/brdc/$@.Z -O $@.Z
	uncompress $@.Z
