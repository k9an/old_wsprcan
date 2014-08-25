/*
 This file is part of k9an-wsprd.
 
 File name: fano.h
 Description:
 
 Copyright 2014, Steven Franke, K9AN
 License: GNU GPL v3
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 */

int fano(unsigned long *metric, unsigned long *cycles,
	unsigned char *data,unsigned char *symbols,
	unsigned int nbits,int mettab[2][256],int delta,
	unsigned long maxcycles);
int encode(unsigned char *symbols,unsigned char *data,unsigned int nbytes);
double gen_met(int mettab[2][256],int amp,double noise,double bias,int scale);

extern unsigned char Partab[];



