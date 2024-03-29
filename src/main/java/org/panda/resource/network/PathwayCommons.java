package org.panda.resource.network;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.biopax.paxtools.pattern.miner.SIFType;
import org.panda.resource.FileServer;
import org.panda.resource.ResourceDirectory;
import org.panda.utility.StringUtil;
import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.Graph;
import org.panda.utility.graph.GraphList;
import org.panda.utility.graph.UndirectedGraph;

import java.io.*;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * Serves Pathway Commons SIF graphs.
 *
 * @author Ozgun Babur
 */
public class PathwayCommons extends FileServer
{
	private static PathwayCommons instance;

	public static PathwayCommons get()
	{
		if (instance == null) instance = new PathwayCommons();
		return instance;
	}

	/**
	 * Gets a phospGraph that is a merge of the desired SIF types.
	 */
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
		return ResourceDirectory.get() + File.separator + "PC/";
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

	/**
	 * Returns the phospGraph for the given SIF type.
	 */
	public Graph getSingleGraph(SIFType type)
	{try{
		String edgeType = type.getTag();
		Graph graph = type.isDirected() ?
			new DirectedGraph("Pathway Commons", edgeType) : new UndirectedGraph("Pathway Commons", edgeType);

		Files.lines(Paths.get(getPrivateDirectory() + type.getTag() + ".txt"))
			.map(line -> line.split("\t")).forEach(token -> {
			if (token.length > 2)
			{
				graph.putRelation(token[0], token[1], token[2]);
			}
			else
			{
				graph.putRelation(token[0], token[1]);
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

	/**
	 * Prints phospGraph statistics.
	 */
	private void printNetworkSizes()
	{
		Set<SIFEnum> types = new HashSet<>(Arrays.asList(SIFEnum.values()));
		types.remove(SIFEnum.NEIGHBOR_OF);
		Graph graph = getGraph(types.toArray(new SIFType[0]));
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
		String base = "https://www.pathwaycommons.org/archives/PC2/v12/";
		String partial = "All.hgnc.txt.gz";

		Scanner sc = new Scanner(new URL(base).openStream());
		while (sc.hasNextLine())
		{
			String filename = StringUtil.fetch(sc.nextLine(), partial, "\"", "\"");
			if (filename != null)
			{
				return new String[]{base + filename};
			}
		}
		throw new RuntimeException("Cannot find the current PC " + partial + " file in " + base);
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
