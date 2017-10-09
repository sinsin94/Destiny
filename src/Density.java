
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.*;

import weka.core.Instances;

public class Density extends Instances {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * The total run time
	 */
	public static long totalRuntime;

	/**
	 * The distance measure.
	 */
	public static final int EUCLIDIAN = 0;

	/**
	 * The distance measure.
	 */
	public static final int MANHATTAN = 1;

	/**
	 * The distance measure.
	 */
	int distanceMeasure;

	/**
	 * Data type.
	 */
	public static final int INFORMATION_SYSTEM = 0;

	/**
	 * Data type.
	 */
	public static final int DECISION_SYSTEM = 1;

	/**
	 * Data type.
	 */
	int dataType;

	/**
	 * How many labels are taught.
	 */
	int numTeach;

	/**
	 * How many labels are predicted.
	 */
	int predicted;

	/**
	 * dc
	 */
	double dc;

	/**
	 * RhoIndex
	 */
	int[] ordrho;

	/**
	 * Rho
	 */
	double[] rho;

	/**
	 * Delta
	 */
	double[] delta;

	/**
	 * Priority
	 */
	double[] priority;

	/**
	 * The node index of centers
	 */
	int[] centers;

	/**
	 * The maximal distance between any pair of points.
	 */
	double maximalDistance;

	/**
	 * The maximal delta
	 */
	double maximalDelta;

	/**
	 * The cluster information. Which cluster do I belong to?
	 */
	int[] clusterIndices;

	/**
	 * The block information.
	 */
	int[][] blockInformation;

	/**
	 * Who is my master?
	 */
	int[] master;

	/**
	 * Predicted labels.
	 */
	int[] predictedLabels;

	/**
	 * Is respective instance already classified? If so, we do not process it
	 * further.
	 */
	boolean[] alreadyClassified;

	/**
	 * The descendant indices to show the importance of instances in a descendant
	 * order.
	 */
	int[] descendantIndices;

	/**
	 * How many instances are taught for each block.
	 */
	int teachInstancesForEachBlock;

	/**
	 ********************************** 
	 * Read from a reader
	 ********************************** 
	 */
	public Density(Reader paraReader) throws IOException, Exception {
		super(paraReader);
		initialize();
	}// Of the first constructor

	/**
	 *************** 
	 * Get the numTeach.
	 **************** 
	 */
	public double getnumTeach() {
		return numTeach;
	}// Of getRi

	/**
	 ********************************** 
	 * Initialize.
	 ********************************** 
	 */
	public void initialize() {
		setDistanceMeasure(EUCLIDIAN);
		setDataType(DECISION_SYSTEM);
		computeMaximalDistance();
	}// Of initialize

	/**
	 ********************************** 
	 * Set Dc from maximal distance.
	 ********************************** 
	 */
	public void setDcFromMaximalDistance(double paraPercentage) {
		dc = maximalDistance * paraPercentage;
		System.out.println("当前dc " + dc);
	}// Of setDcFromMaximalDistance

	/**
	 ********************************** 
	 * Compute distance between instances.
	 ********************************** 
	 */
	public double distance(int paraI, int paraJ, int paraMeasure) {
		distanceMeasure = paraMeasure;
		return distance(paraI, paraJ);
	}// Of distance

	/**
	 ********************************** 
	 * Set the distance measure.
	 ********************************** 
	 */
	public void setDistanceMeasure(int paraMeasure) {
		distanceMeasure = paraMeasure;
	}// Of setDistanceMeasure

	/**
	 ********************************** 
	 * Set the data type.
	 ********************************** 
	 */
	public void setDataType(int paraDataType) {
		dataType = paraDataType;
		setClassIndex(numAttributes() - 1);
	}// Of setDataType

	/**
	 ********************************** 
	 * Compute the maximal distance.
	 ********************************** 
	 */
	public double computeMaximalDistance() {
		maximalDistance = 0;
		double tempDistance;

		for (int i = 0; i < numInstances(); i++) {

			for (int j = 0; j < numInstances(); j++) {
				tempDistance = distanceManhattan(i, j);
				if (maximalDistance < tempDistance) {
					maximalDistance = tempDistance;
				} // Of if
			} // Of for j
		} // Of for i
			// System.out.println("maxdistance are " + maximalDistance);
		System.out.println("最大距离" + maximalDistance);
		return maximalDistance;
	}// Of setDistanceMeasure

	/**
	 ********************************** 
	 * Compute distance between instances.
	 ********************************** 
	 */
	public double distance(int paraI, int paraJ) {
		double tempDistance = 0;
		switch (distanceMeasure) {
		case EUCLIDIAN:
			tempDistance = euclidian(paraI, paraJ);
			break;
		case MANHATTAN:
			tempDistance = manhattan(paraI, paraJ);
			break;
		}// Of switch

		return tempDistance;
	}// Of distance

	/**
	 ********************************** 
	 * Compute distance between instances.
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
	}// Of distance

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

	/**
	 ********************************** 
	 * Manhattan distance.
	 ********************************** 
	 */
	public double manhattan(int paraI, int paraJ) {
		double tempDistance = 0;
		int tempNumAttributes = numAttributes();
		if (dataType == DECISION_SYSTEM) {
			tempNumAttributes--;
		} // Of if
		for (int i = 0; i < tempNumAttributes; i++) {
			tempDistance += Math.abs(instance(paraI).value(i) - instance(paraJ).value(i));
		} // Of for i

		return tempDistance;
	}// Of manhattan

	/**
	 ********************************** 
	 * Manhattan distance of the data set.
	 ********************************** 
	 */
	public double[][] manhattanData() {
		double[][] tempDistance = new double[numInstances()][numInstances()];

		for (int i = 0; i < numInstances(); i++) {

			for (int j = 0; j < numInstances(); j++) {

				tempDistance[i][j] = distanceManhattan(i, j);

			} // of for j
		} // of for i

		return tempDistance;

	}// Of manhattan

	/**
	 ********************************** 
	 * Set dc
	 ********************************** 
	 */
	/**
	 * public void setDc(double paraDc) { dc = paraDc; }// Of setDc
	 * 
	 * /**
	 ********************************** 
	 * Compute rho
	 ********************************** 
	 */
	public void computeRho() {
		rho = new double[numInstances()];

		for (int i = 0; i < numInstances(); i++) {
			for (int j = 0; j < numInstances(); j++) {
				if (distance(i, j) < dc) {
					rho[i]++;
				} // Of if
			} // Of for j
		} // Of for i
	}// Of computeRho

	/**
	 ********************************** 
	 * Compute rho（密度）计算每个实例的密度,该方法很巧妙
	 ********************************** 
	 */
	public void computeRhoWM() {
		rho = new double[numInstances()];

		for (int i = 0; i < numInstances() - 1; i++) {
			for (int j = i + 1; j < numInstances(); j++) {
				if (distance(i, j) < dc) {
					rho[i] = rho[i] + 1;
					rho[j] = rho[j] + 1;
				} // Of if
			} // Of for j
		} // Of for i
	}// Of computeRho

	/**
	 ********************************** 
	 * Compute delta计算距离 找出密度比i大且离i最近的实例作为master[i]=j
	 ********************************** 
	 */
	public void computeDelta() {
		delta = new double[numInstances()];
		master = new int[numInstances()];
		double tempDistance;
		double tempValue;
		int tempMaster;

		for (int i = 0; i < numInstances(); i++) {
			// Initialize. I'm my own master.
			tempMaster = -1;
			tempDistance = Integer.MAX_VALUE;

			for (int j = 0; j < numInstances(); j++) {
				if (rho[j] > rho[i]) {
					tempValue = distance(i, j);
					if (tempValue < tempDistance) {
						tempDistance = tempValue;
						tempMaster = j;
					} // Of if
				} // Of if
			} // Of for j

			delta[i] = tempDistance;
			master[i] = tempMaster;// 如果master[i]==-1,那么i密度最大
									// 如果master[i]是一个非负数，则存在比i密度大的实例，再在比i密度大的实例中考虑距离最小
		} // Of for i

		// Only one real master
		int tempRealMaster = -1;
		for (int i = 0; i < numInstances(); i++) {
			if (master[i] == -1) {// i密度最大
				tempRealMaster = i;
				break;
			} // Of if
		} // Of for i

		for (int i = tempRealMaster + 1; i < numInstances(); i++) {
			System.out.println("master为" + i);
			if (master[i] == -1) {
				master[i] = tempRealMaster;
			} // Of if
		} // Of for i

		// System.out
		// .println("To build the master tree, you need the following array:\r\n"
		// + Arrays.toString(master));
	}// Of computeDelta

	/**
	 ********************************** 
	 * Compute delta
	 ********************************** 
	 */
	public void computeDeltaWM() {
		delta = new double[numInstances()];
		master = new int[numInstances()];
		ordrho = new int[numInstances()];
		ordrho = mergeSortToIndices(rho);// 对密度进行排序
		// System.out.println("ordrho" + Arrays.toString(ordrho));
		delta[ordrho[0]] = -1;
		// 为什么将密度最大的点的master设置为0；密度最大的实例的master不重要，因为并不会用到他的master，它就是老大

		master[ordrho[0]] = 0;
		// 计算每个实例的master，存储的均是坐标，delta数组存储的是实例i与master的距离
		for (int i = 1; i < numInstances(); i++) {// 与比自己密度大的实例进行比较
			delta[ordrho[i]] = maximalDistance;
			for (int j = 0; j <= i - 1; j++) {
				if (distanceManhattan(ordrho[i], ordrho[j]) < delta[ordrho[i]]) {
					delta[ordrho[i]] = distanceManhattan(ordrho[i], ordrho[j]);
					master[ordrho[i]] = ordrho[j];
				} // of if
			} // of for j
		} // of for

		int[] deltaIndex = mergeSortToIndices(delta);

		delta[ordrho[0]] = delta[deltaIndex[0]];// 设置密度最大的点与master的距离为最大
	}// Of computeDelta

	/**
	 ********************************** 
	 * Compute priority. Element with higher priority is more likely to be selected
	 * as a cluster center. Now it is rho * delta. It can also be rho^alpha * delta.
	 ********************************** 
	 */
	public void computePriority() {
		priority = new double[numInstances()];
		for (int i = 0; i < numInstances(); i++) {
			priority[i] = rho[i] * delta[i];
		} // Of for i
	}// Of computePriority

	/**
	 ********************************** 
	 * Compute centers
	 ********************************** 
	 */
	public void computeCenters(int paraNumCenters) {// paraNumCenters = 3;
		centers = new int[paraNumCenters + 1];
		double[] tempWeightedArray = new double[paraNumCenters + 1];// 初始化时，所有的值均为0.0
		double tempValue;
		for (int i = 0; i < numInstances(); i++) {
			tempValue = rho[i] * delta[i];
			for (int j = 0; j < paraNumCenters; j++) {
				if (tempValue > tempWeightedArray[j]) {
					// Move the others
					for (int k = paraNumCenters; k > j; k--) {
						centers[k] = centers[k - 1];
						tempWeightedArray[k] = tempWeightedArray[k - 1];
					} // Of for k
						// Insert here
					centers[j] = i;
					tempWeightedArray[j] = tempValue;
					// Already inserted
					break;
				} // Of if
			} // Of for j
		} // Of for i
	}// Of computeCenters

	/**
	 ********************************** 
	 * Compute centers
	 ********************************** 
	 */
	public void computeCentersWM(int paraNumCenters) {
		centers = new int[paraNumCenters];
		int[] priorityIndex = new int[numInstances()];

		priorityIndex = mergeSortToIndices(priority);
		// System.out.println("this is the test");
		// System.out.println("priorityIndex" + Arrays.toString(priorityIndex));
		for (int i = 0; i < paraNumCenters; i++) {
			centers[i] = priorityIndex[i];
		} // of for i
		System.out.println("centers" + Arrays.toString(centers));
	}// Of computeCenters

	/**
	 ********************************** 
	 * Compute block information
	 ********************************** 
	 */
	public int[][] computeBlockInformation() {
		int tempBlocks = centers.length;
		blockInformation = new int[tempBlocks][];

		for (int i = 0; i < tempBlocks; i++) {
			// Scan to see how many elements
			int tempElements = 0;
			for (int j = 0; j < numInstances(); j++) {
				if (clusterIndices[j] == centers[i]) {// 计算每一簇的数量
					tempElements++;
				} // Of if
			} // Of for k

			// Copy to the list
			blockInformation[i] = new int[tempElements];

			tempElements = 0;
			for (int j = 0; j < numInstances(); j++) {
				if (clusterIndices[j] == centers[i]) {
					blockInformation[i][tempElements] = j;// 记录的是下标
					tempElements++;
				} // Of if
			} // Of for k
		} // Of for i

		System.out.println("blockInformation" + blockInformation.length);

		return blockInformation;
	}// Of computeBlockInformation

	/**
	 ********************************** 
	 * Compute classification accuracy.
	 ********************************** 
	 */
	/**
	 * public double computeClassificationAccuracy() throws Exception { if (dataType
	 * != DECISION_SYSTEM) { throw new Exception( "Classification is inapplicable
	 * for data type: " + dataType); }// Of if
	 * 
	 * int tempCenters = centers.length - 1; int tempClassLabels =
	 * attribute(numAttributes() - 1).numValues(); int[][] tempClusterClassMatrix =
	 * new int[tempCenters][tempClassLabels];
	 * 
	 * int tempClusterCenter, tempClassLabel; for (int i = 0; i < numInstances();
	 * i++) { tempClusterCenter = clusterIndices[i]; tempClassLabel = (int)
	 * instance(i).value(numAttributes() - 1); for (int j = 0; j < tempCenters; j++)
	 * { if (tempClusterCenter == centers[j]) {
	 * tempClusterClassMatrix[j][tempClassLabel]++; }// Of if }// Of for j }// Of
	 * for i
	 * 
	 * // System.out.println(Arrays.deepToString(tempClusterClassMatrix));
	 * 
	 * double tempCorrect = 0; double tempCurrentClusterMaximal; for (int i = 0; i <
	 * tempCenters; i++) { tempCurrentClusterMaximal = 0; for (int j = 0; j <
	 * tempClassLabels; j++) { if (tempCurrentClusterMaximal <
	 * tempClusterClassMatrix[i][j]) { tempCurrentClusterMaximal =
	 * tempClusterClassMatrix[i][j]; }// Of if }// Of for j
	 * 
	 * tempCorrect += tempCurrentClusterMaximal; }// Of for i
	 * 
	 * return tempCorrect / numInstances(); }// Of computeClassificationAccuracy
	 * 
	 * /**
	 ********************************** 
	 * Cluster according to the centers
	 ********************************** 
	 */
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

	/**
	 ********************************** 
	 * Cluster according to the centers
	 * 
	 * @param tempBlocks
	 *            划分的簇的数量
	 ********************************** 
	 */
	public void clusterWithCentersWM(int paraLength) {

		int tempBlocks = paraLength;
		int[] cl = new int[numInstances()];

		clusterIndices = new int[numInstances()];
		for (int i = 0; i < numInstances(); i++) {
			cl[i] = -1;
		} // of for i

		int tempNumber = 0;
		int Ncluster = 0;
		// 初始化，遍历center就好
		// for (int i = 0; i < numInstances(); i++) {
		// if (tempNumber < tempBlocks) {
		// cl[centers[i]] = Ncluster;
		// Ncluster++;
		// tempNumber++;
		// } // of if
		//
		// } // of for i
		for (int i = 0; i < centers.length; i++) {
			cl[centers[i]] = i;
			// Ncluster++;
		}

		// System.out.println("master" + Arrays.toString(master));
		System.out.println("centers" + Arrays.toString(centers));// centers是安装优先级排序的，
		// System.out.println("This is the test 1");
		// System.out.println("cl" + Arrays.toString(cl));
		// System.out.println(ordrho[0]+","+cl[ordrho[0]]);
		// System.out.println(ordrho[1]+","+cl[ordrho[1]]);
		// System.out.println(ordrho[2]+","+cl[ordrho[2]]);
		for (int i = 1; i < numInstances(); i++) {
			// System.out.println(i);
			if (cl[ordrho[i]] == -1) {// 说明不是当前迭代的中心点
				// System.out.println("This is the test 1.1");
				cl[ordrho[i]] = cl[master[ordrho[i]]];// 寻找master是 在比自己密度大的实例中寻找的，好巧妙噢，密度大容易聚在一起
				// System.out.println("当前i是："+i+" ordrho[i]是："+ordrho[i]+"
				// master[ordrho[i]]:"+master[ordrho[i]]+" cl[ordrho[i]]是"+cl[ordrho[i]]);

				// System.out.println("This is the test 1.2");
				// System.out.println("rhoIndexArray"+Arrays.toString(ordrho)+"Now we will
				// label"+i+"labelArray"+Arrays.toString(cl));
				// System.out.println("This is the test 1.2");
			} // of if
		} // of for i
		System.out.println("cl" + Arrays.toString(cl));
		// System.out.println("This is the test 2");
		for (int i = 0; i < numInstances(); i++) {
			// System.out.println("centers[cl[i]];"+centers[cl[i]]);
			clusterIndices[i] = centers[cl[i]];
		}

		// System.out.println("clusterIndices" +
		// Arrays.toString(clusterIndices));
	}// Of clusterWithCenters

	/**
	 ********************************** 
	 * Compute the maximal of a array.
	 ********************************** 
	 */
	public static int getMax(int[] arr) {

		int max = arr[0];

		for (int x = 1; x < arr.length; x++) {
			if (arr[x] > max)
				max = arr[x];

		}
		return max;

	}

	/**
	 ********************************** 
	 * Compute the maximal index of a array.
	 ********************************** 
	 */
	public static int getMaxIndex(double[] arr) {

		double maxNumber = arr[0];
		int maxIndex = 0;
		int tempIndex = 0;

		for (int i = 0; i < arr.length; i++) {
			if (arr[i] > maxNumber) {
				maxNumber = arr[i];
				tempIndex = i;
			}
		}
		maxIndex = tempIndex;

		return maxIndex;

	}

	/**
	 ********************************** 
	 * Cluster based active learning.
	 * 
	 * @param paraBlockStepLength
	 *            The number of blocks added for each round.
	 * @param paraTeachInstances
	 *            How many instances should be taught for each block.
	 ********************************** 
	 */
	/**
	 * public void clusterBasedActiveLearning(int paraBlockStepLength, int
	 * paraTeachInstances) { predictedLabels = new int[numInstances()]; numTeach =
	 * 0; predicted = 0; teachInstancesForEachBlock = paraTeachInstances; boolean[]
	 * tempClassified = new boolean[numInstances()]; boolean[] tempBlockProcessed;
	 * int[][] tempBlockInformation = null; int tempUnprocessedBlocks; int
	 * tempBlocks = paraBlockStepLength; computeRho(); //
	 * System.out.println("clusterBasedActiveLearning test 1"); computeDelta(); //
	 * System.out.println("clusterBasedActiveLearning test 2"); while (true) { //
	 * Step 1. Cluster computeCenters(tempBlocks); clusterWithCenters();
	 * computeBlockInformation(); // System.out.println("clusterBasedActiveLearning
	 * test 2.1");
	 * 
	 * // Step 2. Teach tempUnprocessedBlocks = 0; tempBlockProcessed = new
	 * boolean[tempBlocks]; for (int i = 0; i < tempBlocks; i++) { //
	 * System.out.println("clusterBasedActiveLearning test 2.1.1"); // Which blocks
	 * are processed. tempBlockProcessed[i] = true; for (int j = 0; j <
	 * numInstances(); j++) { if (clusterIndices[j] == centers[i]) { if
	 * (!tempClassified[j]) { tempBlockProcessed[i] = false;
	 * tempUnprocessedBlocks++; break; }// Of if }// Of if }// Of for j //
	 * System.out.println("clusterBasedActiveLearning test 2.1.2"); }// Of for i
	 * 
	 * // All blocks are processed. if (tempUnprocessedBlocks == 0) { break; }// Of
	 * if
	 * 
	 * // Obtain unprocessed block information. tempBlockInformation = new
	 * int[tempUnprocessedBlocks][]; int tempCurrentBlock = 0; //
	 * System.out.println("clusterBasedActiveLearning test 2.2: " // +
	 * tempUnprocessedBlocks); for (int i = 0; i < tempBlocks; i++) { // This block
	 * should not be processed further. if (tempBlockProcessed[i]) { continue; }//
	 * Of if
	 * 
	 * // System.out.println("clusterBasedActiveLearning test 2.2.1"); // Scan to
	 * see how many elements tempBlockInformation[tempCurrentBlock] =
	 * blockInformation[i]; // System.out.println("clusterBasedActiveLearning test
	 * 2.2.2"); tempCurrentBlock++; }// Of for i
	 * 
	 * // System.out.println("clusterBasedActiveLearning test 2.3"); // Teach to
	 * obtain enough labeled instances for (int i = 0; i <
	 * tempBlockInformation.length; i++) { // This block is too small, teach all int
	 * tempUnlabeledInThisBlock = 0; for (int j = 0; j <
	 * tempBlockInformation[i].length; j++) { if
	 * (!tempClassified[tempBlockInformation[i][j]]) { tempUnlabeledInThisBlock++;
	 * }// Of if }// Of for j if (tempUnlabeledInThisBlock <=
	 * teachInstancesForEachBlock) { for (int j = 0; j <
	 * tempBlockInformation[i].length; j++) { if
	 * (!tempClassified[tempBlockInformation[i][j]]) {
	 * tempClassified[tempBlockInformation[i][j]] = true; System.out.println("Teach
	 * " + tempBlockInformation[i][j]); numTeach++; }// Of if }// Of for k }// Of if
	 * // System.out.println("clusterBasedActiveLearning test 2.3.1");
	 * 
	 * int tempTaught = 0; int[] tempRandomIndices =
	 * generateRandomSequence(tempBlockInformation[i].length); for (int j = 0; j <
	 * tempBlockInformation[i].length; j++) { //
	 * System.out.println("clusterBasedActiveLearning test 2.3.1.1"); if
	 * (tempClassified[tempBlockInformation[i][tempRandomIndices[j]]]) { continue;
	 * }// Of if
	 * 
	 * int tempInstanceID = tempBlockInformation[i][tempRandomIndices[j]]; // Set to
	 * the true value. // System.out.println("clusterBasedActiveLearning test
	 * 2.3.1.2"); predictedLabels[tempInstanceID] = (int) instance(
	 * tempInstanceID).classValue(); tempClassified[tempInstanceID] = true;
	 * tempTaught++; numTeach++; System.out.println("Teach " + tempInstanceID + " at
	 * step 2.3"); // System.out.println("clusterBasedActiveLearning test 2.3.1.3");
	 * if (tempTaught >= teachInstancesForEachBlock) { break; }// Of if }// Of for j
	 * }// Of for i // System.out.println("clusterBasedActiveLearning test 2.4");
	 * 
	 * // Step 3. Classify boolean tempPure; int tempCurrentInstance; for (int i =
	 * 0; i < tempBlockInformation.length; i++) { tempPure = true; int tempLabel =
	 * 0;
	 * 
	 * // Is the block pure? boolean tempFirstLabeled = true; for (int j = 0; j <
	 * tempBlockInformation[i].length; j++) { tempCurrentInstance =
	 * tempBlockInformation[i][j]; if (tempClassified[tempCurrentInstance]) { if
	 * (tempFirstLabeled) { tempLabel = predictedLabels[tempCurrentInstance];
	 * tempFirstLabeled = false; } else { if (tempLabel !=
	 * predictedLabels[tempCurrentInstance]) { tempPure = false; break; }// Of if
	 * }// Of if }// Of if }// Of for j
	 * 
	 * // Classify unlabeled instances. if (tempPure) { for (int j = 0; j <
	 * tempBlockInformation[i].length; j++) { tempCurrentInstance =
	 * tempBlockInformation[i][j]; if (!tempClassified[tempCurrentInstance]) {
	 * predictedLabels[tempCurrentInstance] = tempLabel;
	 * tempClassified[tempCurrentInstance] = true; // System.out.println("Classify "
	 * + // tempCurrentInstance); predicted++; }// Of if }// Of if }// Of if }// Of
	 * for i // System.out.println("clusterBasedActiveLearning test 2.5");
	 * 
	 * // Step 4. More blocks. tempBlocks += paraBlockStepLength; }// Of while
	 * 
	 * // System.out.println("clusterBasedActiveLearning finish!");
	 * System.out.println("numTeach = " + numTeach + "; predicted = " + predicted);
	 * System.out.println("Accuracy = " + getPredictionAccuracy()); //
	 * //System.out.println("classified = " + // Arrays.toString(tempClassified));
	 * }// Of clusterBasedActiveLearning
	 * 
	 * /**
	 ********************************** 
	 * Cluster based active learning.
	 * 
	 * @param paraBlockStepLength
	 *            The number of blocks added for each round.7
	 * @param paraTeachInstances
	 *            How many instances should be taught for each block.3
	 * @param paraTeach
	 *            How many instances should be taught.154
	 ********************************** 
	 */
	public void clusterBasedActiveLearningWithSorting(int paraBlockStepLength, int paraTeachInstances, int paraTeach) {
		// Labels for accuracy computation.
		predictedLabels = new int[numInstances()];
		alreadyClassified = new boolean[numInstances()];
		numTeach = 0;
		predicted = 0;
		teachInstancesForEachBlock = paraTeachInstances;
		boolean[] tempBlockProcessed;
		// int[][] tempBlockInformation = null;
		int tempUnprocessedBlocks;
		int tempBlocks = paraBlockStepLength;
		computeRhoWM();
		// System.out.println("clusterBasedActiveLearning test 1");
		computeDeltaWM();
		// System.out.println("clusterBasedActiveLearning test 2");
		computePriority();
		descendantIndices = mergeSortToIndices(priority);
		// System.out
		// .println("The priority of instances in descendant order: \r\n"
		// + Arrays.toString(descendantIndices));
		// while (numTeach + tempBlocks*paraTeachInstances < paraTeach)
		while (true) {
			// Step 2.1 Cluster
			computeCentersWM(tempBlocks);// 根据优先级排序

			clusterWithCentersWM(tempBlocks); // 聚类
			System.out.println("clusterBasedActiveLearning test 2.1");
			// 将cl[]一维数组分成单独的块
			computeBlockInformation();

			// Step 2.2 See which blocks are processed.
			// 当tempUnprocessedBlocks = 0时是一个退出while循环的条件
			tempUnprocessedBlocks = 0;
			tempBlockProcessed = new boolean[tempBlocks];
			for (int i = 0; i < tempBlocks; i++) {
				// System.out.println("clusterBasedActiveLearning test 2.1.1");
				tempBlockProcessed[i] = true;
				for (int j = 0; j < numInstances(); j++) {
					if (clusterIndices[j] == centers[i]) {
						if (!alreadyClassified[j]) {
							tempBlockProcessed[i] = false;
							tempUnprocessedBlocks++;
							break;
						} // Of if
					} // Of if
				} // Of for j
					// System.out.println("clusterBasedActiveLearning test 2.1.2");
			} // Of for i

			// Step 2.3 Finish if all blocks are processed.
			if (tempUnprocessedBlocks == 0) {
				break;
			} // Of if

			// System.out.println("clusterBasedActiveLearning test 2.3");
			// Step 2.4 Check all blocks to teach so that enough labels are
			// obtained
			for (int i = 0; i < blockInformation.length; i++) {// blockInformation是一个二维数组，存聚类信息
				// Step 2.5 No need to teach

				if (tempBlockProcessed[i]) {
					continue;
				} // Of if
					// Step 2.6 This block is too small, teach all.
				if (blockInformation[i].length <= teachInstancesForEachBlock) {
					for (int j = 0; j < blockInformation[i].length; j++) {
						if (!alreadyClassified[blockInformation[i][j]]) {
							if (numTeach > paraTeach) {
								break;
							} // of if
							alreadyClassified[blockInformation[i][j]] = true;
							predictedLabels[blockInformation[i][j]] = (int) instance(blockInformation[i][j])
									.classValue();
							numTeach++;
						} // Of if
					} // Of for k
				} // Of if
					// System.out.println("clusterBasedActiveLearning test 2.3.1");

				// Step 2.7 Organize unknown instances in descendant order
				// according to their priorities.
				int tempTaught = 0;
				int[] tempDescendantIndices = new int[blockInformation[i].length];
				double[] tempPriority = new double[blockInformation[i].length];
				int tempIndex = 0;
				// 将每一块按照密度降序排列
				for (int j = 0; j < numInstances(); j++) {
					if (clusterIndices[descendantIndices[j]] == centers[i]) {// descendantIndices为priority的降序排列
						tempDescendantIndices[tempIndex] = descendantIndices[j];// 存储的依然是坐标，descendantIndices是对实例优先级的排序
						tempPriority[tempIndex] = priority[j];
						tempIndex++;
						// System.out.print(" " + tempIndex);
					} // Of if
				} // Of for j
				for (int j = 0; j < blockInformation[i].length; j++) {
					if (alreadyClassified[tempDescendantIndices[j]]) {
						continue;
					} // Of if
					int tempInstanceID = tempDescendantIndices[j];
					// Set to the true value.
					// System.out.println("clusterBasedActiveLearning test 2.3.1.2");
					if (numTeach > paraTeach) {
						break;
					} // of if
					predictedLabels[tempInstanceID] = (int) instance(tempInstanceID).classValue();
					alreadyClassified[tempInstanceID] = true;
					tempTaught++;
					numTeach++;
					// System.out.print(" " + tempInstanceID);
					// System.out.println("clusterBasedActiveLearning test 2.3.1.3");
					if (tempTaught >= teachInstancesForEachBlock) {
						break;
					} // Of if
				} // Of for j
					// System.out.println();
					// System.out.println("clusterBasedActiveLearningWithSorting test 2.3.2");
					// System.out.println("end, i = " + i);
			} // Of for i
				// System.out.println("clusterBasedActiveLearning test 2.4");

			// Step 3. Classify

			boolean tempPure;
			int tempCurrentInstance;
			for (int i = 0; i < blockInformation.length; i++) {

				if (tempBlockProcessed[i]) {
					continue;
				} // Of if

				tempPure = true;
				int tempLabel = 0;

				// Is the block pure?
				boolean tempFirstLabeled = true;
				for (int j = 0; j < blockInformation[i].length; j++) {
					tempCurrentInstance = blockInformation[i][j];
					if (alreadyClassified[tempCurrentInstance]) {
						if (tempFirstLabeled) {
							tempLabel = predictedLabels[tempCurrentInstance];
							tempFirstLabeled = false;
						} else {
							if (tempLabel != predictedLabels[tempCurrentInstance]) {
								tempPure = false;
								break;
							} // Of if
						} // Of if
					} // Of if
				} // Of for j

				// Classify unlabeled instances.
				if (tempPure) {
					int tempPureNum = 0;
					for (int j = 0; j < blockInformation[i].length; j++) {
						if (alreadyClassified[blockInformation[i][j]]) {
							tempPureNum++;
						} // of if
					} // of for j
					if (tempPureNum > Math.sqrt(blockInformation[i].length)) {
						for (int j = 0; j < blockInformation[i].length; j++) {
							tempCurrentInstance = blockInformation[i][j];
							if (!alreadyClassified[tempCurrentInstance]) {
								predictedLabels[tempCurrentInstance] = tempLabel;
								alreadyClassified[tempCurrentInstance] = true;
								// System.out.print(" " + tempCurrentInstance);
								predicted++;
							} // Of if
						} // Of for j
					} // of if

				} // Of if
			} // Of for i
				// System.out.println("clusterBasedActiveLearning test 2.5");

			// Step 4. More blocks.
			tempBlocks += paraBlockStepLength;
		} // Of while

		// Step 5. imprue, voting

		int max = getMax(predictedLabels);// max = 2;

		computeCenters(max + 1);

		clusterWithCenters();

		computeBlockInformation();

		tempUnprocessedBlocks = 0;
		
		for (int i = 0; i < blockInformation.length; i++) {
			double[] vote = new double[max + 1];
			// System.out.println("test vote");
			for (int j = 0; j < blockInformation[i].length; j++) {
				
				for (int k = 0; k <= max; k++) {
					if (predictedLabels[blockInformation[i][j]] == k) {
						vote[k]++;
					} // of if
				} // of for k
			} // of j
			int voteIndex = 0;
			// System.out.println("blockinformation are "
			// + Arrays.deepToString(blockInformation));

			// System.out.println("vote are" + Arrays.toString(vote));

			voteIndex = getMaxIndex(vote);
			System.out.println("dddddddd"+voteIndex);
			// System.out.println("voteIndex are " + voteIndex);
			for (int j = 0; j < blockInformation[i].length; j++) {
				if (!alreadyClassified[blockInformation[i][j]]) {
					predictedLabels[blockInformation[i][j]] = voteIndex;
					predicted++;
				} // of if

			} // of j

		} // of i

		System.out.println("clusterBasedActiveLearning finish!");
		System.out.println("numTeach = " + numTeach + "; predicted = " + predicted);
		System.out.println("Accuracy = " + getPredictionAccuracy());
		// //System.out.println("classified = " +
		// Arrays.toString(tempClassified));
	}// Of clusterBasedActiveLearningWithSorting

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
	 ******************* 
	 * toString.
	 ******************* 
	 */
	public String toString() {
		String tempString = "\r\nThe center indices are: " + Arrays.toString(centers);
		tempString = "\r\nThe cluster indices are: " + Arrays.toString(clusterIndices);
		return super.toString() + tempString;
	}// Of toString

	/**
	 ******************* 
	 * Density test.
	 ******************* 
	 */
	public static void densityTest() {
		totalRuntime = 0;
		String arffFilename = "/Users/dengsiyu/eclipse/workspace/smaleCode/arff/iris.arff";

		try {
			FileReader fileReader = new FileReader(arffFilename);
			Density tempData = new Density(fileReader);
			fileReader.close();
			tempData.setDistanceMeasure(MANHATTAN);// EUCLIDIAN MANHATTAN
			// System.out.println("This is the data:\r\n" + tempData);

			tempData.setDcFromMaximalDistance(0.2);// 根据两个实例间的最远距离计算衡量距离dc

			tempData.clusterBasedActiveLearningWithSorting(4, 3, 70);

			System.out.println("master are " + Arrays.toString(tempData.master));
			// System.out.println("delta are " +
			// Arrays.toString(tempData.delta));
			// System.out.println("priority are "
			// + Arrays.toString(tempData.priority));

		} catch (Exception ee) {
			// System.out.println("Error occurred while trying to read \'"
			// + arffFilename + "\' in densityTest().\r\n" + ee);
		} // Of try
	}// Of densityTest

	/**
	 ********************************** 
	 * Generate a random sequence of [0, n - 1].
	 * 
	 * @author Hengru Zhang, Revised by Fan Min 2013/12/24
	 * 
	 * @param paraLength
	 *            the length of the sequence
	 * @return an array of non-repeat random numbers in [0, paraLength - 1].
	 ********************************** 
	 */
	public static int[] generateRandomSequence(int paraLength) {
		Random random = new Random();
		// Initialize
		int[] tempResultArray = new int[paraLength];
		for (int i = 0; i < paraLength; i++) {
			tempResultArray[i] = i;
		} // Of for i

		// Swap some elements
		int tempFirstIndex, tempSecondIndex, tempValue;
		for (int i = 0; i < paraLength / 2; i++) {
			tempFirstIndex = random.nextInt(paraLength);
			tempSecondIndex = random.nextInt(paraLength);

			// Really swap elements in these two indices
			tempValue = tempResultArray[tempFirstIndex];
			tempResultArray[tempFirstIndex] = tempResultArray[tempSecondIndex];
			tempResultArray[tempSecondIndex] = tempValue;
		} // Of for i

		return tempResultArray;
	}// Of generateRandomSequence

	/**
	 ********************************** 
	 * Merge sort in descendant order to obtain an index array. The original array
	 * is unchanged.<br>
	 * Examples: input [1.2, 2.3, 0.4, 0.5], output [1, 0, 3, 2].<br>
	 * input [3.1, 5.2, 6.3, 2.1, 4.4], output [2, 1, 4, 0, 3].
	 * 
	 * @author Fan Min 2016/09/09
	 * 
	 * @param paraArray
	 *            the original array
	 * @return The sorted indices.
	 ********************************** 
	 */
	public static int[] mergeSortToIndices(double[] paraArray) {
		int tempLength = paraArray.length;
		int[][] resultMatrix = new int[2][tempLength];
		// Initialize
		int tempIndex = 0;
		for (int i = 0; i < tempLength; i++) {
			resultMatrix[tempIndex][i] = i;
		} // Of for i
			// System.out.println("Initialize, resultMatrix = " +
			// Arrays.deepToString(resultMatrix));

		// Merge
		int tempCurrentLength = 1;
		// The indices for current merged groups.
		int tempFirstStart, tempSecondStart, tempSecondEnd;
		while (tempCurrentLength < tempLength) {
			// System.out.println("tempCurrentLength = " + tempCurrentLength);
			// Divide into a number of groups
			// Here the boundary is adaptive to array length not equal to 2^k.
			// ceil是向上取整函数
			for (int i = 0; i < Math.ceil(tempLength + 0.0 / tempCurrentLength) / 2; i++) {
				// Boundaries of the group
				tempFirstStart = i * tempCurrentLength * 2;
				tempSecondStart = tempFirstStart + tempCurrentLength;// 用于判断是否是最后一块
				tempSecondEnd = tempSecondStart + tempCurrentLength - 1;
				if (tempSecondEnd >= tempLength) {
					tempSecondEnd = tempLength - 1;
				} // Of if
					// System.out.println("tempFirstStart = " + tempFirstStart +
					// ", tempSecondStart = " + tempSecondStart
					// + ", tempSecondEnd = " + tempSecondEnd);

				// Merge this group
				int tempFirstIndex = tempFirstStart;
				int tempSecondIndex = tempSecondStart;
				int tempCurrentIndex = tempFirstStart;
				// System.out.println("Before merge");
				if (tempSecondStart >= tempLength) {
					for (int j = tempFirstIndex; j < tempLength; j++) {
						resultMatrix[(tempIndex + 1) % 2][tempCurrentIndex] = resultMatrix[tempIndex % 2][j];
						tempFirstIndex++;
						tempCurrentIndex++;
					} // Of for j
					break;
				} // Of if

				while ((tempFirstIndex <= tempSecondStart - 1) && (tempSecondIndex <= tempSecondEnd)) {
					if (paraArray[resultMatrix[tempIndex % 2][tempFirstIndex]] >= paraArray[resultMatrix[tempIndex
							% 2][tempSecondIndex]]) {
						// System.out.println("tempIndex + 1) % 2"+tempIndex);
						resultMatrix[(tempIndex + 1) % 2][tempCurrentIndex] = resultMatrix[tempIndex
								% 2][tempFirstIndex];
						// System.out.println("Copy " + resultMatrix[tempIndex %
						// 2][tempFirstIndex] + " from " + tempFirstIndex +
						// " of the first array");
						tempFirstIndex++;
					} else {
						resultMatrix[(tempIndex + 1) % 2][tempCurrentIndex] = resultMatrix[tempIndex
								% 2][tempSecondIndex];
						// System.out.println("Copy " + resultMatrix[tempIndex %
						// 2][tempSecondIndex] + " from " + tempSecondIndex +
						// " of the first array");
						tempSecondIndex++;
					} // Of if
					tempCurrentIndex++;

				} // Of while
					// System.out.println("After compared merge");

				// Remaining part
				// System.out.println("Copying the remaining part");
				for (int j = tempFirstIndex; j < tempSecondStart; j++) {
					resultMatrix[(tempIndex + 1) % 2][tempCurrentIndex] = resultMatrix[tempIndex % 2][j];
					tempCurrentIndex++;
				} // Of for j
				for (int j = tempSecondIndex; j <= tempSecondEnd; j++) {
					resultMatrix[(tempIndex + 1) % 2][tempCurrentIndex] = resultMatrix[tempIndex % 2][j];
					tempCurrentIndex++;
				} // Of for j
					// System.out.println("After copying remaining part");
			} // Of for i
				// System.out.println("Round " + tempIndex + ", resultMatrix = "
				// + Arrays.deepToString(resultMatrix));

			tempCurrentLength *= 2;
			tempIndex++;
		} // Of while
		System.out.println("resultSortedIndices = " + Arrays.toString(resultMatrix[tempIndex % 2]));

		return resultMatrix[tempIndex % 2];
	}// Of mergeSortToIndices

	public static void main(String[] args) {

		long start = System.currentTimeMillis();
		densityTest();

		long end = System.currentTimeMillis();
		System.out.println("���㻨��ʱ��" + (end - start) + "����!");

		// double[] tempArray = {1.2, 2.3, 0.4, 0.5};
		// double[] tempArray = {3.1, 5.2, 6.3, 2.1, 4.4};
		// int[] tempSortedIndices = mergeSortToIndices(tempArray);
		// System.out.println("The indices are: " +
		// Arrays.toString(tempSortedIndices));
	}// Of main
}// Of Density
