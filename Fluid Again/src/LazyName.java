import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;

import javax.swing.JFrame;
import javax.swing.JPanel;

public class LazyName extends JPanel{
	int width = 50, height = 50;
	int ppg = 5;//pixels per gridspace
	double dx = 1; //assuming square gridspaces
	double dt = 0.1;
	int resolution = 100;
	//double viscosity = 0.5;
	double[] pOld, pNew;
	double[] velXOld, velXNew;
	double[] velYOld, velYNew;
	double[] bVelX, bVelY;
	double[] tOld,tNew;
	double[] divergence, pressure;
	boolean[] b;
	
	double time = 0;
	
	double fireConstant = 250;
	
	public static void main(String[] args) {
		JFrame frame = new JFrame("Fluid Redo");
		LazyName fluid = new LazyName();
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.setContentPane(fluid);
		frame.pack();
		frame.setVisible(true);
		for(;;)
		{
			fluid.timeStep();
			fluid.repaint();
			frame.repaint();
			try {
				Thread.sleep(10);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
	}
	
	int I(int x, int y)
	{
//		if(x < 0 || y < 0) 
//			{
//				System.out.println("Out of Bounds: " + x + " " + y);
//				return -1;
//			}
//		if(x >= width || y >= height) 
//			{
//				System.out.println("Out of Bounds: " + x + " " + y);
//				return -1;
//			}
		while(x < 0) x = width - x;
		if(x >= width)  x = x % width;
		while(y < 0) y = height - y;
		if(y >= height) y = y % height;
		return y*width + x;
	}
	
	public LazyName() {
		setPreferredSize(new Dimension(width*ppg,height*ppg));
		init();
	}
	
	public void init()
	{
		pOld = new double[width*height];
		pNew = new double[width*height];
		velXOld = new double[width*height];
		velXNew= new double[width*height];
		velYOld = new double[width*height];
		velYNew = new double[width*height];
		tOld = new double[width*height];
		tNew = new double[width*height];
		b = new boolean[width*height];
		bVelX = new double[width*height];
		bVelY = new double[width*height];
		divergence = new double[width*height];
		pressure = new double[width*height];
		
//		for(int i = 0; i < width*height; i++)
//		{
//			
//		}
		
		for(int x = 0; x < width; x++)
		{
			for(int y = 0; y < height; y++)
			{
				if(x == 0) b[I(x,y)] = true;
				if(y == 0) b[I(x,y)] = true;
				if (x == width-1) b[I(x,y)] = true;
				if(y == height-1) b[I(x,y)] = true;
				
				pOld[I(x,y)] = 10;
				tOld[I(x,y)] = 100;
//				if(x == y)
//				{
//					pOld[I(x,y)] += 30;
//					velXOld[I(x,y)] = 2*x;
//					velYOld[I(x,y)] = -2*y;
//				}
				
			}
		}
	}
	
	public void applyBoundary()
	{
		for(int x = 0; x < width; x++)
		{
			for(int y = 0; y < height; y++)
			{
				velXOld[I(x,y)] = clamp(velXOld[I(x,y)],-10,10);
				velYOld[I(x,y)] = clamp(velYOld[I(x,y)],-10,10);
				pOld[I(x,y)] = clamp(pOld[I(x,y)],0,100);
				//tOld[I(x,y)] = clamp(tOld[I(x,y)],50,150);
			}
		}
	}
	
	public void diffuse(double[] out, double[] in, double diff)
	{
		jacobi_iterative(out,in,1,4);
	}
	
	private void buoyancy(double[] velY) {
		double avgTemp = 0;
		for(int i = 0; i < width*height; i++)
		{
			avgTemp += tOld[i];
		}
		avgTemp /= width*height;
		avgTemp = 1/avgTemp;
		
		if(avgTemp == 0) return; 
		
		double force = 0;
		for(int x = 0; x < width; x++)
		{
			for(int y = 0; y < height; y++)
			{
				if(tOld[I(x,y)] == 0) return; 
				force = avgTemp-(1/tOld[I(x,y)]);
				force *= fireConstant;
				
				velY[I(x,y)] -= force*dt/dx;
			}
		}
	}
	
	public void timeStep()
	{
		time += dt;
		//velXOld[I(25,25)] = 5;
		//velYOld[I(5,5)] = 10;
		
		//applyBoundary();
		
		if(time <= 50)
		{
			velXOld[I(25,25)] = 0.5;
			//velYOld[I(25,25)] = -0.5;
			for(int dx = -2; dx < 2; dx++)
			{
				for(int dy = -2; dy < 2; dy++)
				{
					pOld[I(25+dx,25+dy)] = 50;
				}
			}
			
		}
		
		
		
		advect(pOld,pNew,velXOld,velYOld);
		
		swapP();
		//diffuse(pNew,pOld,1);
		//swapP();
		
		
		advect(tOld,tNew,velXOld,velYOld);
		swapT();
		//jacobi_iterative(tNew,tOld);
		//swapT();
		
		advect(velXOld,velXNew,velXOld,velXOld);
		advect(velYOld,velYNew,velXOld,velXOld);
		
		swapVelX(); swapVelY();
		
		//jacobi_iterative(velXNew,velXOld);
		//jacobi_iterative(velYNew,velYOld);
		
		//buoyancy(velYNew);
		
		project(velXOld,velYOld,divergence,pressure);
		
		//swapVelX(); swapVelY();
		
		
	}
	
    private static double lerp(double s, double e, double t) {
        return s + (e - s) * t;
    }
 
    private static double blerp(double c00, double c10, double c01, double c11, double tx, double ty) {
        return lerp(lerp(c00, c10, tx), lerp(c01, c11, tx), ty);
    }

	public double bilin_interp(double[] grid, double x, double y)
	{
		
		//assuming the coordinates are in gridspace coords
		int x1 = (int) x;
		int y1 = (int) y;
		int x2 = x1 + 1;
		int y2 = y1 + 1;
		
		double dx = x - x1, idx = 1-dx;
		double dy = y - y1, idy = 1-dy;
		
		double q11 = grid[I(x1,y1)];
		double q21 = grid[I(x2,y1)];
		double q12 = grid[I(x1,y2)];
		double q22 = grid[I(x2,y2)];
		
		//double fxy1 = idx*q11 + dx*q21;
		//double fxy2 = idx*q12 + dx*q22;
		
		return blerp(q11,q21,q12,q22,dx,dy);
		
		//return idy*fxy1 + dy*fxy2;
	}
	
	public void advect(double[] in, double[] out, double[] velocityX, double[] velocityY)
	{
		double transform = dt/dx;
		for(int x = 0; x < width; x++)
		{
			for(int y = 0; y < height; y++)
			{
				if(b[I(x,y)]) out[I(x,y)] = 0;
				//TODO transforming from real space to grid space
				//check if this is actually how you do it
				double vx = velocityX[I(x,y)]*transform;
				double vy = velocityY[I(x,y)]*transform;
				
				double px = (x+.5)-vx;
				double py = (y)-vy;
				
				
				
				out[I(x,y)] = bilin_interp(in,px,py);
			}
		}
	}
	
	Color c;
	
	private double clamp(double in, double min, double max)
	{
		if(in < min) return min;
		if(in > max) return max;
		return in;
	}
	
	public void jacobi_iterative(double[] out, double[] in, double a, double c)
	{
		for(int i = 0; i < resolution; i++)
		{
			for(int x = 0; x < width; x++)
			{
				for(int y = 0; y < height; y++)
				{
					double dc = in[I(x,y)];
					double du = out[I(x,y+1)];
					double dr = out[I(x+1,y)];
					double dd = out[I(x,y-1)];
					double dl = out[I(x-1,y)];
					if(b[I(x,y+1)]) du = dc;
					if(b[I(x,y-1)]) dd = dc;
					if(b[I(x+1,y)]) dr = dc;
					if(b[I(x-1,y)]) dl = dc;
					out[I(x,y)] = (du+dr+dd+dl-(a*dc))/c;
				}
			}
		}
	}
	
	public void project(double[] velX, double[] velY, double[] pressure, double[] divergence)
	{
		//calculate divergence and clear the array to hold pressure
		for(int x = 0; x < width; x++)
		{
			for(int y = 0; y < height; y++)
			{
				double du = velY[I(x,y+1)];
				double dr = velX[I(x+1,y)];
				double dd = velY[I(x,y-1)];
				double dl = velX[I(x-1,y)];
				
				if(b[I(x,y+1)]) du = bVelY[I(x,y+1)];
				if(b[I(x,y-1)]) dd = bVelY[I(x,y-1)];
				if(b[I(x+1,y)]) dr = bVelX[I(x+1,y)];
				if(b[I(x-1,y)]) dl = bVelX[I(x-1,y)];
				
				divergence[I(x,y)] = ((dr-dl)+(du-dd))/(2*dx);
				pressure[I(x,y)] = 0;
			}
		}
		
		//solve for pressure
		jacobi_iterative(pressure,divergence,1,4);
		
		for(int x = 0; x < width; x++)
		{
			for(int y = 0; y < height; y++)
			{
				double pu = pressure[I(x,y+1)];
				double pr = pressure[I(x+1,y)];
				double pd = pressure[I(x,y-1)];
				double pl = pressure[I(x-1,y)];
				double pc = pressure[I(x,y)];
				
				double vu = bVelY[I(x,y+1)];
				double vr = bVelX[I(x+1,y)];
				double vd = bVelY[I(x,y-1)];
				double vl = bVelX[I(x-1,y)];
				
				double maskX = 1, maskY = 1;
				double ovX = 0, ovY = 0;
				
				if(b[I(x-1,y)])
				{
					pl = pc;
					ovX = vl;
					maskX = 0;
				}
				if(b[I(x+1,y)])
				{
					pr = pc;
					ovX = vr;
					maskX = 0;
				}
				if(b[I(x,y-1)])
				{
					pd = pc;
					ovY = vd;
					maskY = 0;
				}
				if(b[I(x,y+1)])
				{
					pu = pc;
					ovY = vu;
					maskY = 0;
				}
				
				//subtract the partial derivative of the pressure
				velX[I(x,y)] -= (pr-pl)/(2*dx);
				velY[I(x,y)] -= (pu-pd)/(2*dx);
				
				velX[I(x,y)] *= maskX;
				velY[I(x,y)] *= maskY;
				
				velX[I(x,y)] += ovX;
				velY[I(x,y)] += ovY;
			}
		}
	}
	
	@Override
	protected void paintComponent(Graphics g) {
		super.paintComponent(g);
		for(int x = 0; x < width; x++)
		{
			for(int y = 0; y < height; y++)
			{
				float p = (float) pOld[I(x,y)];
				p = (float) clamp(p,0,100);
				
				float vx = (float) velXOld[I(x,y)];
				//vx = clamp(vx,-5,5);
				
				float vy = (float) velYOld[I(x,y)];
				//vy = clamp(vy,-5,5);
				
				c = new Color(0,0,p/100);
				if(b[I(x,y)]) c = Color.DARK_GRAY;
				g.setColor(c);
				int halfStep = ppg/2;
				g.fillRect(x*ppg, y*ppg, ppg, ppg);
				
				int px = (int) (x*ppg + halfStep);
				int py = (int) (y*ppg + halfStep);
				g.setColor(Color.RED);
				g.drawLine(px, py, px + (int)(vx*ppg/dx), py + (int)(vy*ppg/dx));
				
			}
		}
	}
	
	public void swapP()
	{
		double[] temp = pOld;
		pOld = pNew;
		pNew = temp;
	}
	
	public void swapVelX()
	{
		double[] temp = velXOld;
		velXOld = velXNew;
		velXNew = temp;
	}
	
	public void swapVelY()
	{
		double[] temp = velYOld;
		velYOld = velYNew;
		velYNew = temp;
	}
	
	public void swapT()
	{
		double[] temp = tOld;
		tOld = tNew;
		tNew = temp;
	}
	
}
