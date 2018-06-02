package com.pso;

import java.util.Random;

public class ParticleSwarmOptimization {
	// 權重參數
	private static double w, phi1, phi2;				
	// 最大,小速度限制
	private static double maxVelocity, minVelocity;					
	// 最大,小位置限制
	private static double maxPosition, minPosition;	
	// 粒子數量
	private static int particleCount;					
	// 迭代次數
	private static int iterator;	
	// 粒子
	private static Particle[] particles;
	// 隨機
	private static Random random;
	// 全域最佳解
	private static Particle globalBestParticle;

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		initialParameter();
		initalParticle();
		for(int i = 0; i < iterator; i++) {
			moveParticle();
		}
		System.out.println("done");
	}
	
	static void initialParameter() {
		maxPosition = 100.0;
		minPosition = -100.0;
		maxVelocity = 1.5;
		minVelocity = -1.5;
		w = 1.0;
		phi1 = 0.5;
		phi2 = 0.5;
		particleCount = 100;
		iterator = 5000;
		particles = new Particle[particleCount];
		random = new Random(System.currentTimeMillis());
		globalBestParticle = new Particle();
	}
	
	static void initalParticle() {
		for(int i = 0; i < particleCount; i++) {
			particles[i] = new Particle();
			double position = random.nextDouble() * (maxPosition - minPosition) + minPosition;
			particles[i].position = position;
			particles[i].bestPostion = position;
			double velocity = random.nextDouble() * (maxVelocity - minVelocity) + minVelocity;
			particles[i].velocity = velocity;
			double fitness = fit(particles[i].position);
			particles[i].fintness = fitness;
			particles[i].bestFitness = fitness;
			// 全域最佳設定
			updateGlobalBestFintness(particles[i]);
		}
	}
	
	static void moveParticle() {
		for(int i = 0; i < particleCount; i++) {
			// 更新速度
			double velocity = updateVelocity(particles[i]);
			if(velocity < minVelocity) {
				velocity = minVelocity;
			} else if(velocity > maxVelocity) {
				velocity = maxVelocity;
			}
			// 更新位置
			double position = particles[i].position + velocity;
			if(position < minPosition) {
				position = minPosition;
			} else if(position > maxPosition) {
				position = maxPosition;
			}
			// 更新該粒子狀態
			particles[i].velocity = velocity;
			particles[i].position = position;
			particles[i].fintness = fit(position);
			
			// 更新該粒子目前找過之最佳值
			if(particles[i].fintness < particles[i].bestFitness) {
				particles[i].bestFitness = particles[i].fintness;
				particles[i].bestPostion = particles[i].position;
			}
			
			// 更新全域最佳值
			updateGlobalBestFintness(particles[i]);
		}
	}
	
	static double updateVelocity(Particle particle) {
		//w * v + phi1 * random * (pbest-ppos) + phi2 * random * (gbest-pos);
		 double velocity = w * particle.velocity + 
				 phi1 * random.nextDouble() * (particle.bestPostion - particle.position) + 
				 phi2 * random.nextDouble() * (globalBestParticle.bestPostion - particle.position);
		 return velocity;
	}
	
	static void updateGlobalBestFintness(Particle particle) {
		if(particle.bestFitness < globalBestParticle.bestFitness) {
			globalBestParticle = particle;
//			System.out.println("Now Best Fitness = " + particle.bestFitness);
//			System.out.println("Now Best Position = " + particle.bestPostion);
		}
		System.out.println("Now Best Fitness = " + globalBestParticle.bestFitness);
	}
	
	static double fit(double x) {
		return ackleysFunction(x);
//		return evalAckley(x);
	}
	
	static double ackleysFunction(double x) {		
		double sum1 = 0;
        double sum2 = 0;
        for (int i = 0; i < particleCount; i++) {
            sum1 += x * x;
            sum2 += Math.cos(2 * Math.PI * x);
        }
        double p1 = -20 * Math.exp(-0.2 * Math.sqrt((1.0 / particleCount) * sum1));
        double p2 = Math.exp((1.0 / particleCount) * sum2);
        return p1 - p2 + Math.E + 20;
	}
	
	static double evalAckley(double x) {
		
		double a = 20.0;
		double b = 0.2;
		double c = 2 * Math.PI;
		double firstSum = 0.0;
		double secondSum = 0.0;
		
		for (int i = 0; i < particleCount; i++) {
			firstSum += Math.pow(x, 2);
			secondSum += Math.cos(c * x);
		}
		
		firstSum = -b * Math.sqrt(firstSum / particleCount);
		secondSum = secondSum / particleCount;
		
		return (-a * Math.exp(firstSum)) - secondSum + a + Math.E;
	}
}
