package org.panda.resource;

import org.panda.utility.FileUtil;
import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.Graph;
import org.panda.utility.graph.UndirectedGraph;

import java.io.IOException;
import java.util.*;
import java.util.stream.Stream;

/**
 * This class provides a mapping between HGNC IDs and approved gene symbols.
 *
 * @author Ozgun Babur
 */
public class HGNC extends FileServer
{
	private static HGNC instance;

	private Map<String, String> sym2id;
	private Map<String, String> sym2chr;
	private Map<String, String> id2sym;
	private Map<String, String> old2new;
	private Map<String, String> uniprot2sym;
	private Map<String, String> sym2uniprot;
	private Map<String, Set<String>> families;

	public static synchronized HGNC get()
	{
		if (instance == null) instance = new HGNC();
		return instance;
	}

	public static synchronized void initSingletonWith(Stream<String> resourceStream) throws IOException
	{
		instance = new HGNC()
		{
			@Override
			public synchronized boolean init() throws IOException
			{
				return true;
			}
		};

		instance.load(resourceStream);
	}

	/**
	 * Gets the latest approved official symbol related to the given ID or symbol. If the parameter
	 * is ID, then it should start with "HGNC:".
	 * @param symbolOrID HGNC ID, symbol, or a previous symbol
	 * @return latest symbol
	 */
	public String getSymbol(String symbolOrID)
	{
		if (symbolOrID == null) return null;
		if (id2sym.containsKey(symbolOrID)) return id2sym.get(symbolOrID);
		else if (sym2id.containsKey(symbolOrID)) return symbolOrID;
		symbolOrID = symbolOrID.toUpperCase();
		if (old2new.containsKey(symbolOrID)) return old2new.get(symbolOrID);
		if (uniprot2sym.containsKey(symbolOrID)) return uniprot2sym.get(symbolOrID);
		return null;
	}

	public String getUniProt(String symbolOrID)
	{
		return sym2uniprot.get(getSymbol(symbolOrID));
	}

	public Set<String> getFamily(String name)
	{
		if (!families.containsKey(name)) return Collections.emptySet();
		return families.get(name);
	}

	public Set<String> getAllSymbols()
	{
		return sym2id.keySet();
	}

	public Set<String> getSymbolsOfChromosome(String chr)
	{
		Set<String> set = new HashSet<String>();
		for (String sym : getAllSymbols())
		{
			if (sym2chr.containsKey(sym))
			{
				String c = sym2chr.get(sym);
				String no = c.split("p")[0].split("q")[0];

				String desNo = chr.split("p")[0].split("q")[0];

				if (no.equals(desNo) && c.contains(chr))
				{
					set.add(sym);
				}
			}
		}
		return set;
	}

	public String getChromosomeLoc(String symbol)
	{
		return sym2chr.get(symbol);
	}

	/**
	 * If the given collection
	 */
	public void replaceRecognizedOlds(Collection<String> col, boolean removeUnrecognized)
	{
		Set<String> news = new HashSet<>();
		Set<String> discard = new HashSet<>();
		for (String s : col)
		{
			String n = getSymbol(s);
			if (n == null)
			{
				if (removeUnrecognized) discard.add(s);
			}
			else if (!n.equals(s))
			{
				news.add(n);
				discard.add(s);
			}
		}
		col.removeAll(discard);
		col.addAll(news);
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"hgnc.txt"};
	}

	@Override
	public String[] getDistantURLs()
	{
		return new String[]{"https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym" +
			"&col=gd_prev_sym&col=gd_aliases&col=gd_pub_chrom_map&col=family.name&col=md_prot_id&status=Approved" +
			"&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"};
	}

	@Override
	public boolean load() throws IOException
	{
		return load(getResourceAsStream(getLocalFilenames()[0]));
	}

	public boolean load(Stream<String> resourceStream) throws IOException
	{
		sym2id = new HashMap<>();
		sym2chr = new HashMap<>();
		id2sym = new HashMap<>();
		old2new = new HashMap<>();
		uniprot2sym = new HashMap<>();
		sym2uniprot = new HashMap<>();
		families = new HashMap<>();

		resourceStream.skip(1).forEach(line -> {
			String[] token = line.split("\t");
			String sym = token[1];
			String id = token[0];
			sym2id.put(sym, id);
			id2sym.put(id, sym);

			if (token.length > 2)
			{
				for (String old : token[2].split(","))
				{
					old = old.trim().toUpperCase();
					if (old.isEmpty()) continue;
					old2new.put(old, sym);
				}
			}
			if (token.length > 3)
			{
				for (String synonym : token[3].split(","))
				{
					synonym = synonym.trim().toUpperCase();
					if (synonym.isEmpty()) continue;
					old2new.put(synonym, sym);
				}
			}
			if (token.length > 4)
			{
				sym2chr.put(sym, token[4]);
			}
			if (token.length > 5)
			{
				if (!families.containsKey(token[5]))
					families.put(token[5], new HashSet<>());

				families.get(token[5]).add(sym);
			}
			if (token.length > 6 && !token[6].isEmpty())
			{
				uniprot2sym.put(token[6], sym);
				sym2uniprot.put(sym, token[6]);
			}

		});
		return true;
	}

	public Graph getCompleteClique(boolean directed)
	{
		Graph graph = directed ? new DirectedGraph("HGNC complete clique", "clique-edge") :
			new UndirectedGraph("HGNC complete clique", "clique-edge");

		for (String s1 : sym2id.keySet())
		{
			for (String s2 : sym2id.keySet())
			{
				if (s1.equals(s2)) continue;
				if (!directed && s2.compareTo(s1) < 0) continue;
				graph.putRelation(s1, s2);
				if (directed) graph.putRelation(s2, s1);
			}
		}
		return graph;
	}

	public static void main(String[] args)
	{
//		System.out.println(get().getAllSymbols().size());

//		convertSymbols();
//		convertFile("/home/ozgun/Data/Aslan/platelet-lysine-methylation/140-uniprot.txt");

		String symbol = get().getSymbol("H2AFX");
		System.out.println("symbol = " + symbol);
	}

	private static void convertFile(String file)
	{
		Set<String> notFound = new HashSet<>();

		for (String id : FileUtil.getLinesInSet(file))
		{
			String symbol = get().getSymbol(id);
			if (symbol == null) notFound.add(id);
			else System.out.println(symbol);
		}

		System.out.println("\nNot found:");
		notFound.forEach(System.out::println);
	}

	private static void convertSymbols()
	{
		List<String> list = Arrays.asList(("CaMKIIalpha").split("\n"));

		for (String text : list)
		{
			String sym = get().getSymbol(text);
			System.out.println(text + "\t" + sym);
		}
	}
}
