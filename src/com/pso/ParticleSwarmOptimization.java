package com.pso;

import java.util.Random;

public class ParticleSwarmOptimization {
	// 權重參數
	private static double w, a1, a2;				
	// 最大速度限制
	private static double maxVelocity;					
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
		for(int i = 0; i < iterator; i++ ) {
			moveParticle();
		}
		System.out.println("done");
	}
	
	static void initialParameter() {
		maxPosition = 100.0;
		minPosition = -100.0;
		w = 1.0;
		a1 = 2.0;
		a2 = 2.0;
		particleCount = 20;
		maxVelocity = maxPosition - minPosition; 
		iterator = 100000;
		particles = new Particle[particleCount];
		random = new Random(System.currentTimeMillis());
		globalBestParticle = new Particle();
	}
	
	static void initalParticle() {
		for(int i = 0; i < particleCount; i++) {
			particles[i] = new Particle();
			double position = random.nextFloat() * maxVelocity + minPosition;
			particles[i].position = position;
			particles[i].bestPostion = position;
			particles[i].velocity = random.nextFloat() * maxVelocity;
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
			if(velocity < -maxVelocity) {
				velocity = -maxVelocity;
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
			if(particles[i].fintness > particles[i].bestFitness) {
				particles[i].bestFitness = particles[i].fintness;
				particles[i].bestPostion = particles[i].position;
			}
			
			// 更新全域最佳值
			updateGlobalBestFintness(particles[i]);
		}
	}
	
	static double updateVelocity(Particle particle) {
		//w * v + a1 * random * (pbest-ppos) + a2 * random * (gbest-pos);
		 double velocity = w * particle.velocity + 
				 a1 * random.nextFloat() * (particle.bestPostion - particle.position) + 
				 a2 * random.nextFloat() * (globalBestParticle.bestPostion - - particle.position);
		 return velocity;
	}
	
	static void updateGlobalBestFintness(Particle particle) {
		if(particle.bestFitness > globalBestParticle.bestFitness) {
			globalBestParticle = particle;
			System.out.println("Now Best Fitness = " + particle.bestFitness);
			System.out.println("Now Best Position = " + particle.bestPostion);
		}
	}
	
	static double fit(double x) {
//		return AckleyFunction.ackleysFunction(x);
		return Math.abs(8000.0 + x*(-10000.0+x*(-0.8+x)));
	}
}
