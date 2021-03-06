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
public class ExpressionReader
{
	protected String filename;

	protected Map<String, Map<String, Double>> data;

	protected double LOG2 = Math.log(2);

	protected int idLength;

	public ExpressionReader(String filename) throws FileNotFoundException
	{
		this (filename, null);
	}

	public ExpressionReader(String filename, Set<String> genes) throws FileNotFoundException
	{
		this(filename, genes, 12);
	}

	public ExpressionReader(String filename, Set<String> genes, int idLength) throws FileNotFoundException
	{
		this.filename = filename;
		this.data = new HashMap<>();
		this.idLength = idLength;
		load(genes);
	}

	protected void load(Set<String> genes) throws FileNotFoundException
	{
		Scanner sc = new Scanner(new File(filename));
		String line = sc.nextLine();
		while (line.startsWith("#")) line = sc.nextLine();

		String[] header = line.split("\t");

		// Remove possible quotes
		for (int i = 0; i < header.length; i++)
		{
			if (header[i].startsWith("\"") && header[i].endsWith("\""))
			{
				header[i] = header[i].substring(1, header[i].length() - 1);
			}
		}

		int ss = 0;
		while (!header[ss].startsWith("TCGA")) ss++;

		for (int i = ss; i < header.length; i++)
		{
			if (header[i].length() > idLength) header[i] = header[i].substring(0, idLength);
		}

		// skip second line
		sc.nextLine();

		while (sc.hasNextLine())
		{
			line = sc.nextLine();
			String id = line.substring(0, line.indexOf("|"));
			if (id.startsWith("\"")) id = id.substring(1);

			if (genes != null && !genes.contains(id)) continue;
			if (id.equals("?")) continue;

			String[] token = line.split("\t");

			for (int i = ss; i < header.length; i++)
			{
				if (token[i].equals("NA")) token[i] = "NaN";

				Double val = Double.parseDouble(token[i]);

				if (!data.containsKey(id)) data.put(id, new HashMap<>());
				data.get(id).put(header[i], val);
			}

			if (genes != null && data.size() == genes.size()) break;
		}
	}

	public Set<String> getSamples()
	{
		Set<String> samples = new HashSet<>();
		for (Map<String, Double> map : data.values())
		{
			samples.addAll(map.keySet());
		}
		return samples;
	}

	public Set<String> getGenes()
	{
		return data.keySet();
	}

	public boolean hasGene(String gene)
	{
		return data.containsKey(gene);
	}

	public double[] getGeneAlterationArray(String id, String[] samples)
	{
		if (data.containsKey(id))
		{
			double[] d = new double[samples.length];

			for (int i = 0; i < samples.length; i++)
			{
				if (data.get(id).containsKey(samples[i])) d[i] = Math.log1p(data.get(id).get(samples[i])) / LOG2;
				else d[i] = Double.NaN;
			}
			return d;
		}
		return null;
	}

	public double getGeneAlteration(String id, String sample)
	{
		if (data.containsKey(id) && data.get(id).containsKey(sample))
		{
			return data.get(id).get(sample);
		}
		return Double.NaN;
	}

	public void printStdevHistogram(double binSize)
	{
		Histogram h1 = new Histogram(binSize);
		Histogram h2 = new Histogram(binSize);
		h1.setBorderAtZero(true);
		String[] samples = getSamples().toArray(new String[0]);
		for (String gene : getGenes())
		{
			double[] vals = getGeneAlterationArray(gene, samples);
			double stdev = Summary.stdev(vals);
			h1.count(stdev);
			double max = Summary.max(vals);
			h2.count(max);
		}
		h1.print();
		h2.print();
	}

	public static void main(String[] args) throws FileNotFoundException
	{
		ExpressionReader reader = new ExpressionReader("/home/ozgun/Documents/TCGA/UVM/UVM.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt");
		System.out.println(reader.getSamples().size());
	}
}
