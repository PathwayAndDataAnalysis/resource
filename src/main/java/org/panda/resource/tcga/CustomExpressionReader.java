package org.panda.resource.tcga;

import org.panda.utility.statistics.Histogram;
import org.panda.utility.statistics.Summary;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

/**
 * Reads and serves a TCGA RNA expression file.
 *
 * @author Ozgun Babur
 */
public class CustomExpressionReader extends ExpressionReader
{
	public CustomExpressionReader(String filename) throws FileNotFoundException
	{
		super(filename, null);
	}

	public CustomExpressionReader(String filename, Set<String> genes) throws FileNotFoundException
	{
		super(filename, genes);
	}

	protected void load(Set<String> genes) throws FileNotFoundException
	{
		Scanner sc = new Scanner(new File(filename));
		String line = sc.nextLine();
		while (line.startsWith("#")) line = sc.nextLine();

		String[] header = line.split("\t");

		int ss = 1;

		while (sc.hasNextLine())
		{
			line = sc.nextLine();
			String id = line.substring(0, line.indexOf("\t"));

			if (genes != null && !genes.contains(id)) continue;

			String[] token = line.split("\t");

			for (int i = ss; i < header.length; i++)
			{
				Double val = Double.parseDouble(token[i]);

				if (!data.containsKey(id)) data.put(id, new HashMap<>());
				data.get(id).put(header[i], val);
			}

			if (genes != null && data.size() == genes.size()) break;
		}
	}
}
