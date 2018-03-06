package org.panda.resource;

import org.panda.utility.statistics.Histogram;
import org.panda.utility.statistics.Summary;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * A general reader and server for expression data or other type of data that has genes on the first column, first row
 * is for sample headers, and values presented as a matrix in a delimited text file.
 *
 * @author Ozgun Babur
 */
public class MatrixOfValuesDatasetReader
{
	protected String filename;

	protected Map<String, Map<String, Double>> data;

	protected boolean logTransform;

	/**
	 * Regexp for recognizing the separators of columns.
	 */
	protected String delimiter;

	protected double LOG2 = Math.log(2);

	public MatrixOfValuesDatasetReader(String filename) throws FileNotFoundException
	{
		this.filename = filename;
		this.data = new HashMap<>();
		this.logTransform = false;
		this.delimiter = "\t";
	}

	public void setLogTransform(boolean logTransform)
	{
		this.logTransform = logTransform;
	}

	public void setDelimiter(String delimiter)
	{
		this.delimiter = delimiter;
	}

	public void load(Set<String> genes) throws IOException
	{
		Optional<String> opt = Files.lines(Paths.get(filename)).filter(l -> !l.startsWith("#")).findFirst();
		if (!opt.isPresent()) throw new RuntimeException("Cannot find file header.");

		String[] header = opt.get().split(delimiter);

		Files.lines(Paths.get(filename)).filter(l -> !l.startsWith("#")).skip(1).map(l -> l.split(delimiter))
			.filter(t -> t.length >= header.length).forEach(t ->
			{
				if (genes != null && !genes.contains(t[0])) return;

				for (int i = 1; i < header.length; i++)
				{
					Double val = Double.parseDouble(t[i]);

					if (logTransform) val = Math.log1p(val) / LOG2;

					if (!data.containsKey(t[0])) data.put(t[0], new HashMap<>());
					data.get(t[0]).put(header[i], val);
				}
			});
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
				if (data.get(id).containsKey(samples[i])) d[i] = data.get(id).get(samples[i]);
				else d[i] = Double.NaN;
			}
			return d;
		}
		return null;
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
	}
}
