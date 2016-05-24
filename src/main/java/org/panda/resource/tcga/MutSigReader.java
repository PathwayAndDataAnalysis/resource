package org.panda.resource.tcga;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

/**
 * @author Ozgun Babur
 */
public class MutSigReader
{
	public static Map<String, Double> readPValues(String dir) throws FileNotFoundException
	{
		Map<String, Double> map = new HashMap<>();
		Scanner sc = new Scanner(new File(dir + File.separator + "scores-mutsig.txt"));
		sc.nextLine();
		while (sc.hasNextLine())
		{
			String[] token = sc.nextLine().split("\t");
			String gene = token[1];
			double pval = Double.parseDouble(token[token.length-2]);
			map.put(gene, pval);
		}
		return map;
	}
}
