import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.*;

import weka.core.Instances;

public class DensityMy extends Instances {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * The distance measure.
	 */
	public static final int EUCLIDIAN = 0;

	/**
	 * The distance measure.
	 */
	public static final int MANHATTAN = 1;

	private int distanceMeasure;

	/**
	 * Data type.
	 */
	public static final int INFORMATION_SYSTEM = 0;

	/**
	 * Data type.
	 */
	public static final int DECISION_SYSTEM = 1;

	private int dataType;

	/**
	 * The maximal distance between any pair of points.
	 */
	private double maximalDistance;

	/**
	 * dc
	 */
	private double dc;

	/**
	 * 密度rho[]
	 */

	private double[] rho;

	/**
	 * master
	 */

	private int[] master;

	/**
	 * distancemin
	 */
	private double[] distanceMin;
	/**
	 * 优先级
	 */
	private double[] priority;

	/**
	 * The descendant indices to show the importance of instances in a descendant
	 * order.优先级排序
	 */
	private int[] descendantpriorityIndices;

	/**
	 * centers
	 */
	private int[] centers;

	/**
	 * The cluster information. Which center do I belong to?
	 */
	private int[] clusterIndices;

	private int[][] blockInformation;

	/**
	 * Is respective instance already classified? If so, we do not process it
	 * further.
	 */
	boolean[] alreadyClassified;
	/**
	 * Predicted labels.
	 */
	int[] predictedLabels;

	int numTeach;
	int predicted;

	public DensityMy(Reader paraReader) throws IOException, Exception {
		super(paraReader);
		initialize();
		// TODO Auto-generated constructor stub
	}

	public void initialize() {
		// TODO Auto-generated method stub
		setDistanceMeasure(EUCLIDIAN);
		setDataType(DECISION_SYSTEM);
		computeMaximalDistance();

	}

	private void computeMaximalDistance() {
		// TODO Auto-generated method stub
		maximalDistance = 0;
		double tempDistance;

		for (int i = 0; i < numInstances(); i++) {
			for (int j = 0; j < numInstances(); j++) {
				tempDistance = distanceManhattan(i, j);// 不能用distance，构造函数先调用！
				// System.out.println("tempDistance"+tempDistance);
				if (tempDistance > maximalDistance) {
					maximalDistance = tempDistance;
				} // of if
			} // of for j
		} // of for i
		System.out.println("the maximalDistance: " + maximalDistance);
		// return maximalDistance;
	}// of computeMaximalDistance

	/**
	 ********************************** 
	 * which measure to computer distane.
	 ********************************** 
	 */
	public double distance(int paraI, int paraJ) {
		double tempDistance = 0;

		switch (distanceMeasure) {
		case MANHATTAN:
			tempDistance = distanceManhattan(paraI, paraJ);
			break;

		case EUCLIDIAN:
			tempDistance = euclidian(paraI, paraJ);
			break;
		}
		return tempDistance;

	}// Of distance

	/**
	 ********************************** 
	 * Manhattan distance.
	 ********************************** 
	 */
	public double distanceManhattan(int paraI, int paraJ) {
		double tempDistance = 0;
		int tempNumAttributes = numAttributes();
		if (dataType == DECISION_SYSTEM) {
			tempNumAttributes--;
		} // Of if
		for (int i = 0; i < tempNumAttributes; i++) {
			tempDistance += Math.abs(instance(paraI).value(i) - instance(paraJ).value(i));
		} // Of for i

		return tempDistance;
	}// Of distanceManhattan

	/**
	 ********************************** 
	 * Euclidian distance.
	 ********************************** 
	 */
	public double euclidian(int paraI, int paraJ) {
		double tempDistance = 0;
		double tempValue;
		int tempNumAttributes = numAttributes();
		if (dataType == DECISION_SYSTEM) {
			tempNumAttributes--;
		} // Of if
		for (int i = 0; i < tempNumAttributes; i++) {
			tempValue = (instance(paraI).value(i) - instance(paraJ).value(i));
			tempDistance += tempValue * tempValue;
		} // Of for i

		tempDistance = Math.sqrt(tempDistance);

		return tempDistance;
	}// Of euclidian

	private void setDataType(int decisionSystem) {
		// TODO Auto-generated method stub
		dataType = decisionSystem;
		setClassIndex(numAttributes() - 1);

	}

	public void setDistanceMeasure(int paraMeasure) {
		distanceMeasure = paraMeasure;
		// System.out.println("distanceMeasure"+distanceMeasure);
	}

	private void setDcFromMaximalDistance(double paraPercentage) {
		// TODO Auto-generated method stub
		dc = maximalDistance * paraPercentage;
		System.out.println("the current dc:" + dc);

	}

	private void computerDesity() {
		// TODO Auto-generated method stub

		rho = new double[numInstances()];
		for (int i = 0; i < numInstances() - 1; i++) {

			for (int j = i + 1; j < numInstances(); j++) {
				double tempDistance = distanceManhattan(i, j);
				if (tempDistance < dc) {
					rho[i]++;
					rho[j]++;
				} // of if

			} // of for j

		} // of for i

	}// of computerDesity

	private void computerPriority() {
		// TODO Auto-generated method stub

		priority = new double[numInstances()];
		for (int i = 0; i < numInstances(); i++) {
			priority[i] = rho[i] * distanceMin[i];
		} // of for i
	}// of compuerPriority

	private void computerMaster() {
		// TODO Auto-generated method stub
		master = new int[numInstances()];// 密度大且距离最近,可以先将密度降序处理
		int[] rhoOrder = new int[numInstances()];
		distanceMin = new double[numInstances()];// 存
		rhoOrder = MergeSortToIndex(rho);

		distanceMin[rhoOrder[0]] = -1;
		master[rhoOrder[0]] = 0;

		// 计算每个实例的master，存储的均是坐标，delta数组存储的是实例i与master的距离
		for (int i = 1; i < numInstances(); i++) {// 与比自己密度大的实例进行比较
			distanceMin[rhoOrder[i]] = maximalDistance;

			for (int j = 0; j <= i - 1; j++) {
				double tempRealdistance = distanceManhattan(rhoOrder[i], rhoOrder[j]);
				if (tempRealdistance < distanceMin[rhoOrder[i]]) {
					distanceMin[rhoOrder[i]] = tempRealdistance;
					master[rhoOrder[i]] = rhoOrder[j];
				} // of if
			} // of for j
		} // of for
		int[] deltaIndex = MergeSortToIndex(distanceMin);
		distanceMin[rhoOrder[0]] = distanceMin[deltaIndex[0]];// 设置密度最大的点与master的距离为最大
		System.out.println("the master is :" + Arrays.toString(master));

	}// computerMaster

	private void compuerCenter(int tempBlockNum) {
		// TODO Auto-generated method stub
		centers = new int[tempBlockNum];
		for (int i = 0; i < tempBlockNum; i++) {
			centers[i] = descendantpriorityIndices[i];
		}
		// System.out.println(Arrays.toString(centers));
	}// of compuerCenter
	
	private void compuerCenterdefalt(int paraNumCenters) {
		// TODO Auto-generated method stub
		centers = new int[paraNumCenters+1];//数组长度为4
		
		double[] tempWeightedArray = new double[paraNumCenters + 1];//初始化时，所有的值均为0.0
		double tempValue;
		for (int i = 0; i < numInstances(); i++) {
			tempValue = rho[i] * distanceMin[i];
			for (int j = 0; j < paraNumCenters; j++) {
				if (tempValue > tempWeightedArray[j]) {//大于数组中的优先级，必然要进行插入
					// Move the others
					for (int k = paraNumCenters; k > j; k--) {//在合适位置腾出位置
						centers[k] = centers[k - 1];
						tempWeightedArray[k] = tempWeightedArray[k - 1];
						
					} // Of for k

					// Insert here
					// //System.out.print("Insert at " + j + " ");
					centers[j] = i;
					tempWeightedArray[j] = tempValue;

					// Already inserted
					break;
				} // Of if
			} // Of for j
		} // Of for i
		
		
		
		
	}//compuerCenterdefalt
	
	public void clusterWithCenters() {
		clusterIndices = new int[numInstances()];
		for (int i = 0; i < numInstances(); i++) {
			clusterIndices[i] = -1;
		} // Of for i

		for (int i = 0; i < centers.length - 1; i++) {
			clusterIndices[centers[i]] = centers[i];
		} // Of for i

		// //System.out.println("clusterWithCenters 1");
		int[] tempPathIndices = new int[numInstances()];
		int tempPathLength = 0;
		int tempCurrentIndex;

		for (int i = 0; i < numInstances(); i++) {
			// //System.out.println("clusterWithCenters 2.1");
			// Already processed.
			if (clusterIndices[i] != -1) {
				continue;
			} // Of if

			tempCurrentIndex = i;
			tempPathLength = 0;

			while (clusterIndices[tempCurrentIndex] == -1) {
				tempPathIndices[tempPathLength] = tempCurrentIndex;
				tempCurrentIndex = master[tempCurrentIndex];
				tempPathLength++;
			} // Of while

			// Set master for the path
			for (int j = 0; j < tempPathLength; j++) {
				clusterIndices[tempPathIndices[j]] = clusterIndices[tempCurrentIndex];
			} // Of for j
		} // Of for i
	}// Of clusterWithCenters
	private void computerClusterResult(int paraBlockLength) {
		// TODO Auto-generated method stub

		int[] eachCluster = new int[numInstances()];
		int[] rhoOrder = MergeSortToIndex(rho);

		for (int i = 0; i < numInstances(); i++) {
			eachCluster[i] = -1;
		} // of f

		// 按照centers[]进行划分块
		for (int i = 0; i < centers.length; i++) {
			eachCluster[centers[i]] = i;
		} // of for i

		for (int i = 0; i < numInstances(); i++) {
			if (eachCluster[rhoOrder[i]] == -1) {// 不是中心点
				eachCluster[rhoOrder[i]] = eachCluster[master[rhoOrder[i]]];// eachCluster存储的是类
			}
		} // of for i

		clusterIndices = new int[numInstances()];
		for (int i = 0; i < numInstances(); i++) {

			clusterIndices[i] = centers[eachCluster[i]];
		}
		// System.out.println("clusterIndices"+Arrays.toString(clusterIndices));

	}// computerClusetrResult

	private int[] MergeSortToIndex(double[] paraArray) {

		int tempLength = paraArray.length;
		int[][] resultMatrix = new int[2][tempLength];
		// 初始化，也是为了记录index
		int tempIndex = 0;
		for (int i = 0; i < paraArray.length; i++) {
			resultMatrix[0][i] = i;
		} // of for i

		int tempcurrentLength = 1;
		while (tempcurrentLength < tempLength) {

			for (int i = 0; i < Math.ceil((tempLength + 0.0) / tempcurrentLength / 2.0); i++) {// 定位到哪一块

				int firstStartIndex = i * tempcurrentLength * 2;
				int secondStartIndex = firstStartIndex + tempcurrentLength;
				int secondEndIndex = secondStartIndex + tempcurrentLength - 1;
				if (secondEndIndex >= tempLength) {
					secondEndIndex = tempLength - 1;
				}
				int firstIndex = firstStartIndex;
				int secondIndex = secondStartIndex;
				int tempCurrentIndex = firstStartIndex;

				if (secondStartIndex >= tempLength) {
					for (int j = firstStartIndex; j < tempLength; j++) {
						resultMatrix[(tempIndex + 1) % 2][tempCurrentIndex] = resultMatrix[tempIndex % 2][j];
						tempCurrentIndex++;
					}
					break;
				}
				while ((firstIndex <= secondStartIndex - 1) && (secondIndex <= secondEndIndex)) {// 真正开始做排序的工作

					if (paraArray[resultMatrix[tempIndex % 2][firstIndex]] >= paraArray[resultMatrix[tempIndex
							% 2][secondIndex]]) {

						resultMatrix[(tempIndex + 1) % 2][tempCurrentIndex] = resultMatrix[tempIndex % 2][firstIndex];

						firstIndex++;
					} else {
						resultMatrix[(tempIndex + 1) % 2][tempCurrentIndex] = resultMatrix[tempIndex % 2][secondIndex];

						secondIndex++;
					} // Of if
					tempCurrentIndex++;

				} // Of while

				// Remaining part
				// System.out.println("Copying the remaining part");
				for (int j = firstIndex; j < secondStartIndex; j++) {
					resultMatrix[(tempIndex + 1) % 2][tempCurrentIndex] = resultMatrix[tempIndex % 2][j];
					tempCurrentIndex++;

				} // Of for j
				for (int j = secondIndex; j <= secondEndIndex; j++) {
					resultMatrix[(tempIndex + 1) % 2][tempCurrentIndex] = resultMatrix[tempIndex % 2][j];
					tempCurrentIndex++;
				} // Of for j

			} // of for i

			tempcurrentLength *= 2;
			tempIndex++;

		} // of while
			// System.out.println("resultSortedIndices = " +
			// Arrays.toString(resultMatrix[tempIndex % 2]));
		return resultMatrix[tempIndex % 2];
	}// MergeSort

	private int[][] computerBlockInformation() {
		int tempBlock = centers.length;
		blockInformation = new int[tempBlock][];

		// 计算每一簇的数量

		for (int i = 0; i < tempBlock; i++) {
			int tempEachBlockpnum = 0;
			for (int j = 0; j < numInstances(); j++) {
				if (clusterIndices[j] == centers[i]) {
					tempEachBlockpnum++;
				}
			}
			blockInformation[i] = new int[tempEachBlockpnum];

			int tempnum = 0;
			for (int j = 0; j < numInstances(); j++) {
				// System.out.println("i="+i);
				// System.out.println("tempnum"+tempnum);
				if (clusterIndices[j] == centers[i]) {
					blockInformation[i][tempnum] = j;
					tempnum++;
				} // of if
			} // of for j
				// System.out.println("blockInformation"+blockInformation.length);
		} // of for i

		return blockInformation;
		// TODO Auto-generated method stub

	}

	/**
	 ********************************** 
	 * Compute the maximal of a array.
	 ********************************** 
	 */
	private int getMax(int[] array) {
		// TODO Auto-generated method stub

		int max = array[0];
		for (int i = 0; i < array.length; i++) {
			if (array[i]>max) {
				max = array[i];
			}
		}

		return max;
	}
	
	
	private int getMaxIndex(double[] array) {
		// TODO Auto-generated method stub
		double max = array[0];
		int maxIndex = 0;
		for (int i = 0; i < array.length; i++) {
			if (array[i]>max) {
				max = array[i];
				maxIndex =i;
			}
		}
		
		
		return maxIndex;
	}
	
	/**
	 ******************* 
	 * Get prediction accuracy.
	 ******************* 
	 */
	public double getPredictionAccuracy() {
		double tempInCorrect = 0;
		// System.out.println("Incorrectly classified instances:");
		for (int i = 0; i < numInstances(); i++) {
			if (predictedLabels[i] != (int) instance(i).classValue()) {
				tempInCorrect++;
				System.out.print("" + i + ", ");
			} // Of if
		} // Of for i
		System.out.println();
		System.out.println("This is the incorrect:\r\n" + tempInCorrect);

		return (numInstances() - numTeach - tempInCorrect) / (numInstances() - numTeach);

	}// Of getPredictionAccuracy


	/**
	 ********************************** 
	 * Cluster based active learning.
	 * 
	 * @param paraBlockStepLength
	 *            The number of blocks added for each round.7
	 * @param paraTeachEachInstances
	 *            How many instances should be taught for each block.3
	 * @param paraTeachAll
	 *            How many instances should be taught.154
	 ********************************** 
	 */

	private void testTreePatition(int paraBlockStepLength, int paraTeachEachInstances, int paraAllTeach) {
		// TODO Auto-generated method stub

		int tempBlock = paraBlockStepLength;
		boolean[] tempBlockProcessed;
		predictedLabels = new int[numInstances()];
		alreadyClassified = new boolean[numInstances()];
		numTeach = 0;
		predicted = 0;

		int tempUnprocessedBlocks;
		computerDesity();// 计算实例的密度
		computerMaster();// 计算距离和master
		computerPriority();
		// priority排序
		descendantpriorityIndices = MergeSortToIndex(priority);
		while (true) {
			// step2.1 cluster
			compuerCenter(tempBlock); // 计算中心点
			computerClusterResult(tempBlock);// 按照中心点聚类
			computerBlockInformation();// 计算每一块的信息，存储的是下标
			// step2.2 判断哪些块没有完成
			tempUnprocessedBlocks = 0;// 没有完成的块的数量
			tempBlockProcessed = new boolean[tempBlock];// 初始化为false
			for (int i = 0; i < tempBlock; i++) {
				tempBlockProcessed[i] = true;
				int tempEachBlockNum = blockInformation[i].length;
				for (int j = 0; j < tempEachBlockNum; j++) {
					if (!alreadyClassified[blockInformation[i][j]]) {
						tempBlockProcessed[i] = false;
						// System.out.println(tempUnprocessedBlocks);
						tempUnprocessedBlocks++;
					} // of if
				} // of for j
			} // of for i

			// 所有块都已经被标记
			// System.out.println(tempUnprocessedBlocks);
			if (tempUnprocessedBlocks == 0) {
				break;
			}
			// 存在有些块还没有完成分类，聚焦blockInformation信息
			for (int i = 0; i < blockInformation.length; i++) {
				// 该块的所有实例均已经完成分类
				if (tempBlockProcessed[i]) {
					continue;
				}
				// 该块太小,直接学习
				if (blockInformation[i].length <= paraTeachEachInstances) {
					for (int j = 0; j < blockInformation[i].length; j++) {
						if (!alreadyClassified[blockInformation[i][j]]) {

							if (numTeach >= paraAllTeach) {
								break;
							} // of if
							alreadyClassified[blockInformation[i][j]] = true;
							// 查询标签

							predictedLabels[blockInformation[i][j]] = (int) instance(blockInformation[i][j])
									.classValue();
							numTeach++;
						} // of if

					} // of for j
				} // of if

				// 块比较大，学习优先级高的一些实例，所以对优先级进行排序
				int tempTeach = 0;
				int[] tempDescendantIndices = new int[blockInformation[i].length];// 用于存储具体的优先级
				// tempDescendantIndices =MergeSortToIndex(blockInformation[i]);
				double[] tempPriorities = new double[blockInformation[i].length];

				int tempIndex = 0;
				// 对每一块的优先级进行排序
				for (int j = 0; j < numInstances(); j++) {
					if (clusterIndices[descendantpriorityIndices[j]] == centers[i]) {
						tempDescendantIndices[tempIndex] = descendantpriorityIndices[j];
						tempPriorities[tempIndex] = priority[j];

						tempIndex++;
					} // of if
				} // of for j

				// 对优先级高的实例进行学习，先判断是否学习过
				for (int j = 0; j < blockInformation[i].length; j++) {
					if (alreadyClassified[tempDescendantIndices[j]]) {
						continue;
					} // of if

					if (numTeach >= paraAllTeach) {
						break;
					} // of if
					predictedLabels[tempDescendantIndices[j]] = (int) instance(tempDescendantIndices[j]).classValue();
					alreadyClassified[tempDescendantIndices[j]] = true;
					// predicted++;
					numTeach++;
					tempTeach++;
					if (tempTeach >= paraTeachEachInstances) {
						break;
					}
				} // of for j

			} // of for i
				// 每一块学习完毕后进行分类预测
			for (int i = 0; i < blockInformation.length; i++) {
				// 该块的所有实例均已经完成分类
				if (tempBlockProcessed[i]) {
					continue;
				}
				// 所学习到的是否纯？
				boolean tempPure = true;
				int templabel = 0;
				boolean currentfirtstLabel = true;
				for (int j = 0; j < blockInformation[i].length; j++) {

					if (alreadyClassified[blockInformation[i][j]]) {
						if (currentfirtstLabel) {
							templabel = predictedLabels[blockInformation[i][j]];
							currentfirtstLabel = false;
						} else {
							if (templabel != predictedLabels[blockInformation[i][j]]) {
								tempPure = false;
								break;
							} // of if
						} // if else

					} // of if

				} // of for j

				// 纯的话学习
				if (tempPure) {
					for (int j = 0; j < blockInformation[i].length; j++) {
						// 将剩余的分为同一类

						if (!alreadyClassified[blockInformation[i][j]]) {
							predictedLabels[blockInformation[i][j]] = templabel;
							alreadyClassified[blockInformation[i][j]] = true;
							predicted++;
						} // of if

					} // of for j
				} // of if

			} // of for i

			tempBlock += paraBlockStepLength;
		} // of while

		// impure voting
		int max = getMax(predictedLabels);
		//System.out.println(max);
		compuerCenterdefalt(max+1);
		clusterWithCenters();
		computerBlockInformation();
		tempUnprocessedBlocks = 0;
		for (int i = 0;i < blockInformation.length; i++) {
			double[] vote = new double[max+1];
			
			for (int j = 0; j < blockInformation[i].length; j++) {
				
				
				for (int k = 0; k < vote.length; k++) {
					if (predictedLabels[blockInformation[i][j]]==k) {						
						vote[k]++;					
					}//of if 
				}//of for k
				
			}//of for j
			
			
			int maxindex = getMaxIndex(vote);
			for (int j = 0; j <blockInformation[i].length; j++) {
				if (!alreadyClassified[blockInformation[i][j]]) {
					predictedLabels[blockInformation[i][j]]=maxindex;
					alreadyClassified[blockInformation[i][j]] = true;
					predicted++;
					
				}
				
			}
			
			
		}//of for i
		
		
		
		System.out.println("clusterBasedActiveLearning finish!");
		System.out.println("numTeach = " + numTeach + "; predicted = " + predicted);
		System.out.println("Accuracy = " + getPredictionAccuracy());

	}// of testTreePatition

	

	

	private static void densityTest() {
		// TODO Auto-generated method stub
		String fileName = "/Users/dengsiyu/eclipse/workspace/smaleCode/arff/iris.arff";

		try {
			FileReader paraReader = new FileReader(fileName);
			DensityMy tempData = new DensityMy(paraReader);

			paraReader.close();

			tempData.setDistanceMeasure(MANHATTAN);// EUCLIDIAN MANHATTAN
			tempData.setDcFromMaximalDistance(0.2);
			tempData.testTreePatition(4, 3, 70);

		} catch (Exception ee) {
			// TODO Auto-generated catch block
			ee.printStackTrace();
		}

	}

	public static void main(String[] args) {

		long start = System.currentTimeMillis();
		densityTest();

		long end = System.currentTimeMillis();
		// System.out.println("锟斤拷锟姐花锟斤拷时锟斤拷" + (end - start) + "锟斤拷锟斤拷!");

		// double[] tempArray = {1.2, 2.3, 0.4, 0.5};
		// double[] tempArray = {3.1, 5.2, 6.3, 2.1, 4.4};
		// int[] tempSortedIndices = mergeSortToIndices(tempArray);
		// System.out.println("The indices are: " +
		// Arrays.toString(tempSortedIndices));
	}// Of main

}
