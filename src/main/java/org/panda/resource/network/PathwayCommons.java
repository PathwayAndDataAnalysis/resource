package org.panda.resource.network;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.biopax.paxtools.pattern.miner.SIFType;
import org.panda.resource.FileServer;
import org.panda.utility.StringUtil;
import org.panda.utility.graph.Graph;
import org.panda.utility.graph.GraphList;

import java.io.*;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class PathwayCommons extends FileServer
{
	public Graph getGraph(SIFType... types)
	{
		if (fileExists(types))
		{
			if (types.length == 1) return getSingleGraph(types[0]);
			else if (types.length > 1)
			{
				GraphList graph = new GraphList("Pathway Commons");

				for (SIFType type : types)
				{
					graph.addGraph(getSingleGraph(type));
				}

				return graph;
			}
		}
		return null;
	}

	protected String getPrivateDirectory()
	{
		return getResourceFilesLocation() + File.separator + "PC/";
	}

	protected boolean fileExists(SIFType... types)
	{
		for (SIFType type : types)
		{
			if (!Files.exists(Paths.get(getPrivateDirectory() + type.getTag() + ".txt")))
				return false;
		}
		return true;
	}

	public Graph getSingleGraph(SIFType type)
	{try{
		String edgeType = type.getTag();
		Graph graph = new Graph("Pathway Commons", edgeType);

		Files.lines(Paths.get(getPrivateDirectory() + type.getTag() + ".txt"))
			.map(line -> line.split("\t")).forEach(token -> {
			if (token.length > 2)
			{
				graph.putRelation(token[0], token[1], token[2], type.isDirected());
			}
			else
			{
				graph.putRelation(token[0], token[1], type.isDirected());
			}
		});

		return graph;
	}
	catch (IOException e){throw new RuntimeException(e);}}


	public boolean processTheDownloadedFiles()
	{try{
		File file = new File(locateInBase(getLocalFilenames()[0]));
		Scanner sc = new Scanner(file);
		sc.nextLine(); // skip header

		Map<String, Writer> writers = new HashMap<>();

		Files.createDirectories(Paths.get(getPrivateDirectory()));

		while (sc.hasNextLine())
		{
			String line = sc.nextLine();
			if (line.isEmpty()) break;

			String[] token = line.split("\t");

			if (token.length > 2)
			{
				if (!writers.containsKey(token[1])) writers.put(token[1],
					new BufferedWriter(new FileWriter(getPrivateDirectory() + token[1] + ".txt")));

				writers.get(token[1]).write(token[0] + "\t" + token[2]);

				if (token.length > 6)
				{
					writers.get(token[1]).write("\t" + token[6]);
				}

				writers.get(token[1]).write("\n");
			}
		}

		for (Writer writer : writers.values())
		{
			writer.close();
		}

		return file.delete();
	}
	catch (IOException e){throw new RuntimeException(e);}}

	public static void main(String[] args)
	{
		PathwayCommons pc = new PathwayCommons();
//		printDataOverlaps();
		pc.printNetworkSizes();
//		printMostConnected(getGraph(SIFEnum.CONTROLS_STATE_CHANGE_OF), 500);
	}

	private void printNetworkSizes()
	{
		Graph graph = getGraph(SIFEnum.values());
		graph.printStats();
	}

	private void printMostConnected(Graph graph, int limit)
	{
		List<String> genes = new ArrayList<String>(graph.getSymbols());
		final Map<String, Integer> degree = new HashMap<String, Integer>();
		for (String gene : genes)
		{
			degree.put(gene, graph.getDegree(gene));
		}
		Collections.sort(genes, (o1, o2) -> degree.get(o2).compareTo(degree.get(o1)));

		int i = 0;
		for (String gene : genes)
		{
			i++;
			System.out.println(gene + "\t" + degree.get(gene));
			if (i == limit) break;
		}
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"PC.sif"};
	}

	@Override
	public String[] getDistantURLs() { try
	{
		String base = "http://www.pathwaycommons.org/archives/PC2/current/";
		String partial = "All.EXTENDED_BINARY_SIF.hgnc.txt.gz";

		Scanner sc = new Scanner(new URL(base).openStream());
		while (sc.hasNextLine())
		{
			String filename = StringUtil.fetch(sc.nextLine(), partial, "\"", "\"");
			if (filename != null)
			{
				return new String[]{base + filename};
			}
		}
		throw new RuntimeException("Cannot find the current PC All.EXTENDED_BINARY_SIF.hgnc.txt.gz file in " + base);
	}
	catch (Exception e){throw new RuntimeException(e);}}

	@Override
	public boolean load() throws IOException
	{
		return true;
	}

	@Override
	public boolean localResourceExists()
	{
		return Files.exists(Paths.get(getPrivateDirectory()));
	}
}
