package com.pso;

public class AckleyFunction {

	static double ackleysFunction(double x) {
		double p1 = -20 * Math.exp(-0.2 * Math.sqrt(0.5 * ((x * x))));
		double p2 = Math.exp(0.5 * (Math.cos(2 * Math.PI * x)));
		return p1 - p2 + Math.E + 20;
	}
}
