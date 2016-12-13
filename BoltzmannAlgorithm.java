/*
	Author: Darren Schlager
*/

import java.util.Scanner;
import java.io.FileReader;
import java.io.IOException;;

public class BoltzmannAlgorithm 
{

	static double[] t = {-0.1,-0.2,0.7};	// thresholds
	static double[][] w = { 				// weights
						{0,-0.5,0.4}, 
						{-0.5,0,0.5}, 
						{0.4,0.5,0} 
	 				};
	
	public static void main(String[] args) 
	{
		initialize(); // set the class variables using the contents of a file; comment out this line to use the default values

		System.out.println("Thresholds\n");
		printMatrix(t, 1);
		System.out.println("\nWeights\n");
		printMatrix(w, 1);
		
		double[][] s = generateStateMatrix(t.length);
		System.out.println("\nStates\n");
		printMatrix(s, 0);
		
		double[][] wMs = matrixMultiply(w, s);
		System.out.println("\nW * S\n");
		printMatrix(wMs, 1);
		
		double[][] wMsSt = matrixSubtract(wMs, t);
		System.out.println("\nW * S - T\n");
		printMatrix(wMsSt, 1);
		
		double[][] hl = matrixHardLimit(wMsSt);
		System.out.println("\nHL[ W * S - T ]\n");
		printMatrix(hl, 0);
		
		double[][] transitions = calculateTransitionMatrix(hl);
		System.out.println("\nTransitions\n");
		printMatrix(transitions, 0);
		
		double[][][] transitionMatrices = generateTransitionMatrices(s, hl);
		for(int i=0; i<transitionMatrices.length; i++) {
			System.out.println();
			for(int j=0; j<s.length; j++)
			{
				System.out.printf("%4.0f",s[j][i]);
			}
			System.out.println("\n |–");
			printMatrix(transitionMatrices[i], 0);		}
		
		double[] energies = calculateEnergies(s, w, t);
		System.out.println("\nEnergies\n");
		for(int i=0; i<s[0].length; i++)
		{
			for(int j=0; j<s.length; j++)
			{
				System.out.printf("%5.0f",s[j][i]);
			}
			System.out.print("   =   ");
			double num = energies[i]*10;
			num = Math.round(num);
			num /= 10;
			System.out.printf("%4.1f\n",num);
		}
		
		double[][][] probabilities = calculateProbabilities(wMsSt);
		System.out.println("\nProbabilities\n");
		for(int i=0; i<probabilities.length; i++)
		{
			for(int j=0; j<s.length; j++)
			{
				System.out.printf("%3.0f",s[j][i]);
			}
			System.out.println("      P[0]     P[1] \n|–");
			for(int j=0; j<probabilities[0].length; j++)
			{
				for(int k=0; k<transitionMatrices[0][0].length; k++)
				{
					System.out.printf("%3.0f", transitionMatrices[i][j][k]);
				}
				System.out.print("  ");
				for(int k=0; k<probabilities[0][0].length; k++)
				{
					System.out.printf("%9.4f", probabilities[i][j][k]);
				}
				System.out.println();
			}
			System.out.println();
		}
		
		double[][] probabilityTransitionMatrix = generateProbabilityTransitionMatrix(probabilities, transitionMatrices, s);
		System.out.print("\nProbability Transition Matrix\n\n");
		for(int i=0; i<probabilityTransitionMatrix.length; i++)
		{
			System.out.printf("%8d",i);
		}
		System.out.println("\n  +–");
		for(int i=0; i<probabilityTransitionMatrix.length; i++)
		{
			System.out.printf("%-3d",i);
			for(int j=0; j<probabilityTransitionMatrix.length; j++)
			{
				if(probabilityTransitionMatrix[i][j]!=0)
				{
					System.out.printf("%7.4f ", probabilityTransitionMatrix[i][j]);
				}
				else 
				{
					System.out.print("   ––   ");
				}
			}
			System.out.println();
		}
		
		double[][] steadyStateVector = matrixMultiply(probabilityTransitionMatrix, probabilityTransitionMatrix);
		double epsilon = 0.000000000001;
		int numIterations = 9999;
		for(int i=0; i<numIterations; i++)
		{
			steadyStateVector = matrixMultiply(steadyStateVector, probabilityTransitionMatrix);
			for(int k=0; k<steadyStateVector.length; k++)
			{
				boolean foundSteadyStateVector = true;
				for(int j=1; j<steadyStateVector.length; j++)
				{
					if(Math.abs(steadyStateVector[0][k]-steadyStateVector[j][k])>epsilon)
					{
						foundSteadyStateVector = false;
						break;
					}
				}
				if(k == steadyStateVector.length-1 && foundSteadyStateVector)
				{
					numIterations = i;
					
				}
			}
		}
		System.out.print("\nSteady State Vector\n\n   ");
		printMatrix(steadyStateVector[0], 4);
		System.out.println("\n   iteration count: " + (numIterations+1));
		
	}
	
	/*
		sets the class variables using the contents of a file
	*/
	public static void initialize()
	{
		// used to retrieve keyboard input from the user
		Scanner keyboard = new Scanner(System.in);
		
		// print description of expeteted file contents
		System.out.println("BoltzmannAlgorithm\n");
		System.out.println("Format your file as follows:");
		System.out.println("===============================================");
		System.out.println(" <number nodes>");
		System.out.println(" <threshold 1> <threshold 2> ... <threshold n>");
		System.out.println(" <weight 1,1> <weight 1,2> ... <weight 1,n>");
		System.out.println(" <weight 2,1> <weight 2,2> ... <weight 2,n>");
		System.out.println(" .");
		System.out.println(" .");
		System.out.println(" .");
		System.out.println(" <weight n,1> <weight n,2> ... <weight n,n>");
		System.out.println("===============================================");
		
		Scanner file;
		boolean fileProcessedSuccessfully = false;
		do {
			
			// prompt the user for the file path
			System.out.print("file path: ");
			String path = keyboard.nextLine();
			System.out.println();
			
			try 
			{
				// valid file
				
				// open the file
				file = new Scanner(new FileReader(path));
				
				//get the number of input nodes
				int numNodes = 0;
				try 
				{
					numNodes = file.nextInt();

					// set the class variable to contain the threshold for each node
					t = new double[numNodes];
					try 
					{
						for(int i=0; i<numNodes; i++)
						{
							t[i] = file.nextDouble();
						}
						
						try 
						{
							// set the class variable to contain the weights of each connection
							w = new double[numNodes][numNodes];
							for(int i=0; i<numNodes; i++)
							{
								for(int j=0; j<numNodes; j++)
								{
									w[i][j] = file.nextDouble();
								}
							}
							
							fileProcessedSuccessfully = true;
							
						} catch (Exception e) // w
						{
							System.out.println("That file is not formatted correctly. Try Again.");
						}
					}
					catch (Exception e) // t
					{
						System.out.println("That file is not formatted correctly. Try Again.");
					}
				}
				catch (Exception e) // numNodes
				{
					System.out.println("That file is not formatted correctly. Try Again.");
				}
				
				//close the file
				file.close();
			} 
			catch (IOException e) // file
			{
				// invalid file
				System.out.println("That file does not exist. Try Again.");
			}
			
		} while (!fileProcessedSuccessfully);
	}
	
	/*
		generates a matrix that contains a binary number represented by each column, from 0 to n-1
	*/
	public static double[][] generateStateMatrix(int num) 
	{
		
		int numStates = 1;
		for(int i=0; i<num; i++) 
		{
			numStates *= 2;
		}
		
		double[][] result = new double[num][numStates];
		int alternateEvery = 1;
		int currentAlternateEvery = 1;
		
		for(int i=num-1; i>=0; i--) 
		{
			
			int currentNum = 0;
			
			for(int j=0; j<numStates; j++)
			{
				result[i][j] = currentNum;
				
				if(--currentAlternateEvery == 0) 
				{
					
					currentAlternateEvery = alternateEvery;
					
					if(currentNum==1) 
					{
						currentNum = 0;
					}
					else 
					{
						currentNum = 1;
					}
				}
			}
			
			currentAlternateEvery = alternateEvery*=2;
		}
		
		return result;
	}
	
	/*
		generates a matrix that contains the list of states in row 1 and which state they transition to in row 2
	*/
	public static double[][] calculateTransitionMatrix(double[][] hl) 
	{
		double[][] result = new double[2][hl[0].length];
		
		for(int i=0; i<result[0].length; i++) {
			result[0][i] = i;
		}
		
		double multiplier = Math.pow(2, hl.length-1);
		for(int i=0; i<hl.length; i++)
		{
			for(int j=0; j<hl[0].length; j++)
			{
				if(hl[i][j]==1)
				{
					result[1][j] += multiplier;
				} 
			}
			multiplier /= 2;
		}
		
		return result;
	}
	
	/*
		generates list of transition matrices, one for each state
			- each matrix is nXn where n is the number of rows in the state matrix (s)
			- each matrix is initialized to contain the the corresponding column in the state matrix (s) repeated in each row; 
			  then the diagonal is replaced by the corresponding column from the hard limit matrix (hl)
			
	*/
	public static double[][][] generateTransitionMatrices(double[][] s, double[][] hl) 
	{
		double[][][] results = new double[s[0].length][][];
		
		for(int i=0; i<results.length; i++)
		{
			double[][] result = new double[s.length][s.length];
			
			for(int j=0; j<s.length; j++)
			{
				for(int k=0; k<s.length; k++)
				{
					result[j][k] = s[k][i];
				}
			}
			
			for(int j=0; j<s.length; j++)
			{
				result[j][j] = hl[j][i];
			}
			
			results[i] = result;
		}
		
		return results;
	}
	
	/*
		generates a array that contains the energies for each state
	*/
	public static double[] calculateEnergies(double[][] s, double[][] w, double[] t) 
	{
		double[] result = new double[s[0].length];
		
		for(int i=0; i<result.length; i++)
		{
			double[] v = new double[s.length];
			for(int j=0; j<s.length; j++)
			{
				v[j] = s[j][i];
			}
			result[i] = calculateEnergy(v, w, matrixTranspose(t));
		}
		
		return result;
	}
	
	/*
		computes the energy of a state
	*/
	public static double calculateEnergy(double[] v, double[][] w, double[][] t)
	{
		double[] result1 = matrixMultiply(v, w);
		double[] result2 = matrixMultiply(result1, matrixTranspose(v));
		for(int i=0; i<result2.length; i++)
		{
			result2[i] *= -0.5;
		}
		double[] result3 = matrixMultiply(v, t);
		double[] result4 = matrixAdd(result2, result3);
		
		return result4[0];
	}
	
	/*
		generates list of probability matrices, one for each state
			- each matrix contains the probability of transitioning to a 0 and to a 1 for each transition in the transition matrix for that state
	*/
	public static double[][][] calculateProbabilities(double[][] wMsST)
	{
		double[][][] result = new double[wMsST[0].length][wMsST.length][2];
		for(int i=0; i<result.length; i++)
		{
			for(int j=0; j<result[0].length; j++)
			{
				result[i][j][1] = 1 / ( 1 + Math.pow(Math.E, (-1*wMsST[j][i])/0.5) );
				result[i][j][0] = 1 - result[i][j][1];
				
				result[i][j][1] /= wMsST.length;
				result[i][j][0] /= wMsST.length;
			}
		}
		return result;
	}
	
	/*
		uses a list of transitions and probability matrices for each state to generate a matrix that contains the probabilities of transitioning from one state to another
	*/
	public static double[][] generateProbabilityTransitionMatrix(double[][][] probabilities, double[][][] transitionMatrices, double[][] s) 
	{
		double[][] result = new double[s[0].length][s[0].length];
		for(int i=0; i<result.length; i++)
		{
			for(int j=0; j<result[0].length; j++)
			{
				result[i][j] = 0.0;
			}
		}
		
		for(int i=0; i<result.length; i++)
		{
			double probabilityOfNotTransitioning = 1.0;
			for(int j=0; j<s.length; j++)
			{
				if(s[j][i] != transitionMatrices[i][j][j])
				{
					if(transitionMatrices[i][j][j]==0) 
					{
						result[i][binaryToDecimal(transitionMatrices[i][j])] = probabilities[i][j][0];
						probabilityOfNotTransitioning -= probabilities[i][j][0];
					}
					else 
					{
						result[i][binaryToDecimal(transitionMatrices[i][j])] = probabilities[i][j][1];
						probabilityOfNotTransitioning -= probabilities[i][j][1];
					}
				}
				else
				{
					double[] state = new double[transitionMatrices[i][j].length];
					System.arraycopy(transitionMatrices[i][j], 0, state, 0, transitionMatrices[i][j].length);
					if(transitionMatrices[i][j][j]==0) 
					{
						state[j] = 1;
						result[i][binaryToDecimal(state)] = probabilities[i][j][1];
						probabilityOfNotTransitioning -= probabilities[i][j][1];
					}
					else 
					{
						state[j] = 0;
						result[i][binaryToDecimal(state)] = probabilities[i][j][0];
						probabilityOfNotTransitioning -= probabilities[i][j][0];
					}
				}
			}
			result[i][i] = probabilityOfNotTransitioning;
		}
		
		return result;
	}
	
	/*
		converts an array containg a binary number to a decimal integer
	*/
	public static int binaryToDecimal(double[] arr)
	{
		int multiplier = 1;
		for(int i=0; i<arr.length-1; i++)
		{
			multiplier *= 2;
		}
		
		int result = 0;
		for(int i=0; i<arr.length; i++)
		{
			if(arr[i]==1)
			{
				result += multiplier;
			}
			multiplier /= 2;
		}
		
		return result;
	}
	
	/*
		matrix hard limit (Nxn)
			- takes the hard limit of each value: (>=0) = 1 , (<0) = 0 
	*/
	public static double[][] matrixHardLimit(double[][] arr) 
	{
		double[][] result = new double[arr.length][arr[0].length];
		
		for(int i=0; i<arr.length; i++)
		{
			for(int j=0; j<arr[0].length; j++)
			{
				if(arr[i][j]<0)
				{
					result[i][j] = 0;
				}
				else 
				{
					result[i][j] = 1;	
				}
			}
		}
		
		return result;
	}
	
	/*
		the row in a Nx1 matrix becomes the column in a 1xN matrix	
	*/
	public static double[][] matrixTranspose(double[] arr)
	{
		double[][] result = new double[arr.length][1];
		for(int i=0; i<arr.length; i++)
		{
			result[i][0] = arr[i];
		}
		return result;
	}
	
	/*
		matrix addition (1xN)
			- perform subtraction on each element with the same index
			- the matrices must have the same dimensions 
	*/
	public static double[] matrixAdd(double[] arr1, double arr2[]) 
	{
		if(arr1.length==arr2.length)
		{
			double[] result = new double[arr1.length];
			for(int i=0; i<arr1.length; i++) 
			{
				result[i] = arr1[i]+arr2[i];
			}
			return result;
		}
		else
		{
			return null;	
		}
	}
	
	/*
		matrix subtraction (MxN - Mx1)
			- subtract the row in the second matrix from each column in the first matrix
			- the number of rows in the first matrix must match the number of columns in the second matrix
	*/
	public static double[][] matrixSubtract(double[][] arr1, double arr2[]) 
	{
		if(arr1.length==arr2.length)
		{
			double[][] result = new double[arr1.length][arr1[0].length];
			for(int j=0; j<arr1[0].length; j++) 
			{
				for(int i=0; i<arr2.length; i++) {
					result[i][j] = arr1[i][j]-arr2[i];
				}
			}
			return result;
		}
		else
		{
			return null;	
		}
	}
	
	/*
		matrix multiplication (1xM, MxP)
			- for each column in the second matrix
				+ multiply each element in the row of the first matrix by the corresponding element in the column; sum the results
			- the number of columns in the first matrix MUST equal the number of rows in the second matrix
				 
	*/
	public static double[] matrixMultiply(double[] arr1, double arr2[][]) 
	{
		if(arr1.length==arr2.length) 
		{
			int numCol = arr2[0].length;
			double[] result = new double[numCol];
			
			for(int i=0; i<numCol; i++)
			{
				double sum = 0;
				for(int j=0; j<arr1.length; j++) 
					sum += arr1[j]*arr2[j][i];
				result[i] = sum;
			}
			
			return result;
		}
		else 
		{
			return null;
		}
	}
	
	/*
		matrix multiplication (NxM, MxP)
			- for each row in the first matrix, repeat the following for each column in the second matrix
				+ multiply each element in the row by the corresponding element in the column; sum the results
			- the number of columns in the first matrix MUST equal the number of rows in the second matrix
				 
	*/
	public static double[][] matrixMultiply(double[][] arr1, double arr2[][]) 
	{
		if(arr1[0].length==arr2.length) 
		{
			
			int numRow = arr1.length;
			int numCol = arr2[0].length;
			double[][] result = new double[numRow][numCol];
			
			for(int i=0; i<numRow; i++) 
			{
				for(int j=0; j<numCol; j++)
				{
					double sum = 0;
					for(int k=0; k<arr1[i].length; k++)
						sum += arr1[i][k]*arr2[k][j];
					result[i][j] = sum;
				}
			}
			
			return result;
		}
		else 
		{
			return null;
		}
	}
	
	/*
		print the contents of a NxN matrix	
	*/
	public static void printMatrix(double[][] matrix, int places) 
	{
		String str = "%";
		str += (4+places)+".";
		str += "" + places;
		str += "f";
		for(int i=0; i<matrix.length; i++) 
		{
			for (int j=0; j<matrix[i].length; j++) 
			{
				System.out.printf(str, matrix[i][j]);
			}
			System.out.println();
		}
	}
	
	/*
		print the contents of a Nx1 matrix
	*/
	public static void printMatrix(double[] matrix, int places)
	{
		String str = "%";
		str += (4+places)+".";
		str += "" + places;
		str += "f";
		for (int i=0; i<matrix.length; i++) 
		{
			System.out.printf(str, matrix[i]);
		}
		System.out.println();
	}
	
}