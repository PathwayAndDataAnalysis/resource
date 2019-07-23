package org.panda.resource.autismdatasets;

import org.panda.resource.FileServer;
import org.panda.resource.HG37;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.Progress;
import org.panda.utility.statistics.Histogram;
import org.panda.utility.statistics.KernelDensityPlot;
import org.panda.utility.statistics.KolmogorovSmirnov;
import org.panda.utility.statistics.UniformityChecker;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Provides genomic alterations in denovo-db.
 *
 * @author Ozgun Babur
 */
public class DenovoDB extends FileServer
{
	private static DenovoDB instance;

	private List<Entry> data;

	private Map<String, String> patientToControlMap;
	private Map<String, String> controlToPatientMap;
	private Set<String> yuen2017T3SampleNames;
	private Set<String> yuen2017LowQualitySampleNames;
	private Set<String> yuenAnCommonGenes;
	private Set<String> yuenTurnerCommonGenes;
	private Set<String> yuenTurnerAnCommonGenes;
	private Set<String> turnerAnCommonGenes;
	private Set<String> turnerSamples;

	public static final DataFilter SELECT_MUT = e -> !e.gene.equals("NA") && !e.gene.isEmpty();
//		&& (!(e.proteinVariant.isEmpty() || e.proteinVariant.equals("NA")) || e.functionClass.startsWith("splice"));
//		&& !(!(e.proteinVariant.isEmpty() || e.proteinVariant.equals("NA")) || e.functionClass.startsWith("splice"))
//		&& !e.functionClass.equals("intron");

//	private static final Set<String> GOOD_AUTISM_STUDIES = new HashSet<>(Arrays.asList("DeRubeis2014", "Hashimoto2015", "Iossifov", "Krumm", "Turner_2017", "Yuen2016", "Yuen2017"));
//	private static final Set<String> GOOD_AUTISM_STUDIES = new HashSet<>(Arrays.asList("DeRubeis2014", "Hashimoto2015", "Iossifov", "Krumm"));
	private static final Set<String> GOOD_AUTISM_STUDIES = new HashSet<>(Arrays.asList( "An2018", "Yuen2017"));
	public static final DataFilter SELECT_AUTISM_STUDY = e -> GOOD_AUTISM_STUDIES.contains(e.studyName) &&
		(!e.studyName.equals("Yuen2017") || get().yuen2017T3SampleNames.contains(e.sampleID));// &&
//		!get().getYuen2017LowQualitySampleNames().contains(e.sampleID);


	public static synchronized DenovoDB get()
	{
		if (instance == null) instance = new DenovoDB();
		return instance;
	}

	public Set<String> getAllGenes()
	{
		return data.stream().map(d -> d.gene).filter(g -> !g.equals("NA") && !g.isEmpty()).collect(Collectors.toSet());
	}

	public List<Entry> getData()
	{
		return data;
	}

	public Stream<Entry> getDataStream(DataFilter filter)
	{
		if (filter != null) return data.stream().filter(filter::select);
		else return data.stream();
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"denovo-db.non-ssc-samples.variants.tsv", "denovo-db.ssc-samples.variants.tsv",
			"Yuen2017-T3.csv", "Joon.txt"};
	}

	@Override
	public String[] getDistantURLs()
	{
		return new String[]{"http://denovo-db.gs.washington.edu/denovo-db.non-ssc-samples.variants.tsv.gz",
			"http://denovo-db.gs.washington.edu/denovo-db.ssc-samples.variants.tsv.gz",
			GITHUB_REPO_BASE + "Yuen2017-T3.csv", GITHUB_REPO_BASE + "Joon.txt"};
	}

	@Override
	public boolean load() throws IOException
	{
		data = new ArrayList<>(900000);

		getResourceAsStream(getLocalFilenames()[0]).skip(2).forEach(l -> data.add(new Entry(l)));
		getResourceAsStream(getLocalFilenames()[1]).skip(2).forEach(l -> data.add(new Entry(l)));
		getResourceAsStream(getLocalFilenames()[3]).skip(2).forEach(l -> data.add(new Entry(l, true)));

		yuen2017T3SampleNames = FileUtil.getTermsInTabDelimitedColumn(locateInBase(getLocalFilenames()[2]), 0, 2);

		DataFilter mssng_has_gene_filter = e -> e.studyName.equals("Yuen2017") && yuen2017T3SampleNames.contains(e.sampleID) && !e.gene.equals("NA") && !e.gene.isEmpty();
		yuenAnCommonGenes = CollectionUtil.getIntersection(
			getDataStream(mssng_has_gene_filter).map(e -> e.gene).collect(Collectors.toSet()),
			getDataStream(DataFilterEnum.AN_HAS_GENE).map(e -> e.gene).collect(Collectors.toSet()));

		turnerAnCommonGenes = CollectionUtil.getIntersection(
			getDataStream(DataFilterEnum.TURNER_HAS_GENE).map(e -> e.gene).collect(Collectors.toSet()),
			getDataStream(DataFilterEnum.AN_HAS_GENE).map(e -> e.gene).collect(Collectors.toSet()));

		yuenTurnerCommonGenes = CollectionUtil.getIntersection(
			getDataStream(mssng_has_gene_filter).map(e -> e.gene).collect(Collectors.toSet()),
			getDataStream(DataFilterEnum.TURNER_HAS_GENE).map(e -> e.gene).collect(Collectors.toSet()));

		yuenTurnerAnCommonGenes = CollectionUtil.getIntersection(
			getDataStream(mssng_has_gene_filter).map(e -> e.gene).collect(Collectors.toSet()),
			getDataStream(DataFilterEnum.TURNER_HAS_GENE).map(e -> e.gene).collect(Collectors.toSet()),
			getDataStream(DataFilterEnum.AN_HAS_GENE).map(e -> e.gene).collect(Collectors.toSet()));

		turnerSamples = getDataStream(DataFilterEnum.TURNER).map(e -> e.sampleID).collect(Collectors.toSet());

		populatePairedControls();
		return true;
	}

	public void plot()
	{
	}

	private void populatePairedControls()
	{
		patientToControlMap = new HashMap<>();
		controlToPatientMap = new HashMap<>();
		Set<String> pSet = new HashSet<>();
		Set<String> sSet = new HashSet<>();

		data.forEach(e ->
		{
			if (e.sampleID.endsWith(".p1"))
			{
				String key = e.sampleID.substring(0, e.sampleID.lastIndexOf("."));
				pSet.add(key);
			}
			else if (e.sampleID.endsWith(".s1"))
			{
				String key = e.sampleID.substring(0, e.sampleID.lastIndexOf("."));
				sSet.add(key);
			}
		});

		Set<String> common = CollectionUtil.getIntersection(pSet, sSet);
		for (String key : common)
		{
			patientToControlMap.put(key + ".p1", key + ".s1");
			controlToPatientMap.put(key + ".s1", key + ".p1");
		}
	}

	public Set<String> getControlledPatients()
	{
		return patientToControlMap.keySet();
	}

	public String getControlOfPatient(String sample)
	{
		if (patientToControlMap.containsKey(sample))
		{
			return patientToControlMap.get(sample);
		}
		return null;
	}

	public String getPatientOfControl(String sample)
	{
		if (controlToPatientMap.containsKey(sample))
		{
			return controlToPatientMap.get(sample);
		}
		return null;
	}

	public Set<String> getControlledSamples()
	{
		Set<String> set = new HashSet<>(patientToControlMap.keySet());
		set.addAll(controlToPatientMap.keySet());
		return set;
	}

	public Set<String> getYuen2017LowQualitySampleNames()
	{
		if (yuen2017LowQualitySampleNames == null)
		{
			Map<String, List<Entry>> map = get().getDataStream(e -> e.studyName.equals("Yuen2017") &&
				e.primaryPhenotype.equals("autism") && get().yuen2017T3SampleNames.contains(e.sampleID) &&
				SELECT_MUT.select(e))
				.collect(Collectors.groupingBy(e -> e.sampleID));

			Map<String, Long> counts = map.keySet().stream().collect(Collectors.toMap(Function.identity(),
				id -> map.get(id).stream().map(e -> e.gene).distinct().count()));

			yuen2017LowQualitySampleNames = counts.keySet().stream().filter(sample -> counts.get(sample) < 18)
				.collect(Collectors.toSet());
		}
		return yuen2017LowQualitySampleNames;
	}

	public class Entry
	{
		public String sampleID;
		public String studyName;
		public String pubmedID;
		public int numProbands;
		public int numControls;
		public String sequenceType;
		public String primaryPhenotype;
		public String validation;
		public String chr;
		public int position;
		public String variant;
		public String rsID;
		public String dbsnpBuild;
		public String ancestralAllele;
		public int genome1000Count;
		public double exacFreq;
		public double espAaFreq;
		public double espEaFreq;
		public String transcript;
		public int codingDnaSize;
		public String gene;
		public String functionClass;
		public String cDnaVariant;
		public String proteinVariant;
		public String exonIntron;
		public double polyPhenHDiv;
		public double polyPhenHVar;
		public double siftScore;
		public double caddScore;
		public double lofScore;
		public double lrtScore;

		public Entry(String line)
		{
			initFromDBLine(line);
		}

		private void initFromDBLine(String line)
		{
			String[] t = line.split("\t");
			sampleID = t[0];
			studyName = t[1];
			pubmedID = t[2];
			numProbands = Integer.valueOf(t[3]);
			numControls = Integer.parseInt(t[4]);
			sequenceType = t[5];
			primaryPhenotype = t[6];
			validation = t[7];
			chr = t[8];
			position = Integer.valueOf(t[9]);
			variant = t[10];
			rsID = t[11];
			dbsnpBuild = t[12];
			ancestralAllele = t[13];
			genome1000Count = Integer.valueOf(t[14]);
			exacFreq = Double.valueOf(t[15]);
			espAaFreq = Double.valueOf(t[16]);
			espEaFreq = Double.valueOf(t[17]);
			transcript = t[18];
			codingDnaSize = Integer.valueOf(t[19]);
			gene = t[20];
			functionClass = t[21];
			cDnaVariant = t[22];
			proteinVariant = t[23];
			exonIntron = t[24];
			polyPhenHDiv = Double.valueOf(t[25]);
			polyPhenHVar = Double.valueOf(t[26]);
			siftScore = Double.valueOf(t[27]);
			caddScore = Double.valueOf(t[28]);
			lofScore = Double.valueOf(t[29]);
			lrtScore = Double.valueOf(t[30]);
		}

		private final static String JOON_STUDY_NAME = "An2018";
		private final static String JOON_PMID = "30545852";
		private final static int JOON_PROBANDS = 1902;
		private final static String JOON_SEQ_TYPE = "genome";
		private final static String JOON_AUTISM = "autism";
		private final static String JOON_CONTROL = "control";
		private final static String JOON_VALIDATION = "unknown";

		public Entry(String line, boolean fromJoon)
		{
			if (!fromJoon)
			{
				initFromDBLine(line);
				return;
			}

			String[] t = line.split("\t");
			sampleID = t[7];
			this.studyName = JOON_STUDY_NAME;
			pubmedID = JOON_PMID;
			numProbands = JOON_PROBANDS;
			numControls = JOON_PROBANDS;
			sequenceType = JOON_SEQ_TYPE;
			primaryPhenotype = t[6].equals("case") ? JOON_AUTISM : JOON_CONTROL;
			validation = JOON_VALIDATION;
			chr = t[0].substring(3);
			position = Integer.valueOf(t[1]);
			variant = t[2] + ">" + t[3];
			transcript =  t[10].equals(".") ? "" : t[10];
			gene =  t[9].equals(".") ? "" : t[9];
			functionClass =  t[8];
			proteinVariant =  "";
		}

		public String key()
		{
			return studyName + sampleID + chr + position;
		}
	}

	public interface DataFilter
	{
		boolean select(Entry e);
	}

	public enum DataFilterEnum implements DataFilter
	{
		HAS_GENE(e -> !e.gene.equals("NA") && !e.gene.isEmpty()),
		INTRON(e -> HAS_GENE.select(e) && e.functionClass.contains("intron")),
		HAS_GENE_BUT_NOT_INTRON(e -> HAS_GENE.select(e) && !e.functionClass.contains("intron")),

		YUEN(e -> e.studyName.equals("Yuen2017") && DenovoDB.get().yuen2017T3SampleNames.contains(e.sampleID)),
		TURNER(e -> e.studyName.equals("Turner_2017")),
		TURNER_REDO(e -> e.studyName.equals("An2018") && DenovoDB.get().turnerSamples.contains(e.sampleID)),
		AN(e -> e.studyName.equals("An2018")),
		AN_DIF(e -> e.studyName.equals("An2018") && !DenovoDB.get().turnerSamples.contains(e.sampleID)),

		AUTISM(e -> e.primaryPhenotype.equals("autism")),
		CONTROL(e -> e.primaryPhenotype.equals("control")),

		YUEN_AUTISM(e -> HAS_GENE.select(e) && YUEN.select(e) && AUTISM.select(e)),
		TURNER_AUTISM(e -> HAS_GENE.select(e) && TURNER.select(e) && AUTISM.select(e)),
		TURNER_REDO_AUTISM(e -> HAS_GENE.select(e) && TURNER_REDO.select(e) && AUTISM.select(e)),
		AN_AUTISM(e -> HAS_GENE.select(e) && AN.select(e) && AUTISM.select(e)),
		AN_CONTROL(e -> HAS_GENE.select(e) && AN.select(e) && CONTROL.select(e)),
		AN_DIF_AUTISM(e -> HAS_GENE.select(e) && AN_DIF.select(e) && AUTISM.select(e)),

		YUEN_HAS_GENE(e -> HAS_GENE.select(e) && YUEN.select(e)),
		TURNER_HAS_GENE(e -> HAS_GENE.select(e) && TURNER.select(e)),
		AN_HAS_GENE(e -> HAS_GENE.select(e) && AN.select(e)),

		YUEN_TURNER_CONTROL(e -> HAS_GENE.select(e) && (YUEN.select(e) || TURNER.select(e)) && CONTROL.select(e)),
		YUEN_AUTISM_TURNER_CONTROL(e -> HAS_GENE.select(e) && ((YUEN.select(e) && AUTISM.select(e))  || (TURNER.select(e) && CONTROL.select(e)))),

		YUEN_AN_AUTISM(e -> HAS_GENE.select(e) && (YUEN.select(e) || AN.select(e)) && AUTISM.select(e)),
		YUEN_AN(e -> HAS_GENE.select(e) && (YUEN.select(e) || AN.select(e))),
		YUEN_TURNER_AUTISM(e -> HAS_GENE.select(e) && (YUEN.select(e) || TURNER.select(e)) && AUTISM.select(e)),
		YUEN_TURNER_REDO_AUTISM(e -> HAS_GENE.select(e) && (YUEN.select(e) || TURNER_REDO.select(e)) && AUTISM.select(e)),
		YUEN_TURNER_AN_AUTISM(e -> HAS_GENE.select(e) && (YUEN.select(e) || TURNER.select(e) || AN_DIF.select(e)) && AUTISM.select(e)),
		YUEN_TURNER_AN_AUTISM_INTRON(e -> INTRON.select(e) && YUEN_TURNER_AN_AUTISM.select(e)),
		YUEN_TURNER_AN_AUTISM_NOT_INTRON(e -> HAS_GENE_BUT_NOT_INTRON.select(e) && YUEN_TURNER_AN_AUTISM.select(e)),
		YUEN_AN_AUTISM_WITH_COMMON_GENES(e -> YUEN_AN_AUTISM.select(e) && DenovoDB.get().yuenAnCommonGenes.contains(e.gene)),
		YUEN_TURNER_AN_AUTISM_WITH_COMMON_GENES(e -> YUEN_TURNER_AN_AUTISM.select(e) && DenovoDB.get().yuenTurnerAnCommonGenes.contains(e.gene)),
		YUEN_TURNER_AUTISM_WITH_COMMON_GENES(e -> YUEN_TURNER_AUTISM.select(e) && DenovoDB.get().yuenTurnerCommonGenes.contains(e.gene)),
		YUEN_AN_CONTROL(e -> HAS_GENE.select(e) && (YUEN.select(e) || AN.select(e)) && CONTROL.select(e)),
		YUEN_AUTISM_AN_CONTROL(e -> HAS_GENE.select(e) && ((YUEN.select(e) && AUTISM.select(e))  || (AN.select(e) && CONTROL.select(e)))),
		TURNER_AN_AUTISM(e -> HAS_GENE.select(e) && (AN_DIF.select(e) || TURNER.select(e)) && AUTISM.select(e)),
		TURNER_AN_AUTISM_WITH_COMMON_GENES(e -> TURNER_AN_AUTISM.select(e) && DenovoDB.get().turnerAnCommonGenes.contains(e.gene)),

		YUEN_TURNER_AUTISM_INTRON(e -> INTRON.select(e) && (YUEN.select(e) || TURNER.select(e)) && AUTISM.select(e)),
		YUEN_TURNER_CONTROL_INTRON(e -> INTRON.select(e) && (YUEN.select(e) || TURNER.select(e)) && CONTROL.select(e)),

		YUEN_AN_AUTISM_INTRON(e -> INTRON.select(e) && (YUEN.select(e) || AN.select(e)) && AUTISM.select(e)),
		YUEN_AN_CONTROL_INTRON(e -> INTRON.select(e) && (YUEN.select(e) || AN.select(e)) && CONTROL.select(e)),

		YUEN_TURNER_AUTISM_NOT_INTRON(e -> HAS_GENE_BUT_NOT_INTRON.select(e) && (YUEN.select(e) || TURNER.select(e)) && AUTISM.select(e)),
		YUEN_TURNER_CONTROL_NOT_INTRON(e -> HAS_GENE_BUT_NOT_INTRON.select(e) && (YUEN.select(e) || TURNER.select(e)) && CONTROL.select(e)),

		YUEN_AN_AUTISM_NOT_INTRON(e -> HAS_GENE_BUT_NOT_INTRON.select(e) && (YUEN.select(e) || AN.select(e)) && AUTISM.select(e)),
		YUEN_AN_CONTROL_NOT_INTRON(e -> HAS_GENE_BUT_NOT_INTRON.select(e) && (YUEN.select(e) || AN.select(e)) && CONTROL.select(e)),
		;

		DataFilter filter;

		DataFilterEnum(DataFilter filter)
		{
			this.filter = filter;
		}

		@Override
		public boolean select(Entry e)
		{
			return filter.select(e);
		}

		public static DataFilterEnum get(String name)
		{
			name = name.toUpperCase().replaceAll("-", "_");
			return valueOf(name);
		}
	}

	public static void main(String[] args) throws IOException
	{
//		printStudiesAndPhenotypeFrequencies();
//		printMutationDensityPerStudy();
//		printSampleOverlaps();
//		printMutationOverlapInCommonSamples("Iossifov", "Krumm");
//		printRowsPerSamplePerStudy();
//		printGenesPerSamplePerStudy();
//		printMutationClassCounts();
//		recordMutationNonUniformity();
//		printMutationUniformityAroundAGene();
//		plotHotspotsForAGene();

		printStudiesContent();

//		temp();
	}

	private static void printStudiesAndPhenotypeFrequencies()
	{
		List<String> phenotypes = get().getDataStream(e -> true).map(e -> e.functionClass).distinct().sorted()
			.collect(Collectors.toList());

		List<String> studies = get().getDataStream(e -> true).map(e -> e.studyName).distinct().sorted()
			.collect(Collectors.toList());

		for (String phenotype : phenotypes)
		{
			System.out.print("\t" + phenotype);
		}

		for (String study : studies)
		{
			System.out.print("\n" + study);

			for (String phenotype : phenotypes)
			{
				long count = get().getDataStream(e -> e.studyName.equals(study) && e.functionClass.equals(phenotype))
//					.map(e -> e.sampleID)
					.filter(SELECT_MUT::select)
					.distinct().count();
				System.out.print("\t" + count);
			}
		}
	}

	private static void printMutationDensityPerStudy()
	{
		List<String> studies = get().getDataStream(e -> e.primaryPhenotype.equals("autism")).map(e -> e.studyName)
			.distinct().sorted().collect(Collectors.toList());

		for (String study : studies)
		{
			long sampleCnt = get().getDataStream(e -> e.studyName.equals(study)).map(e -> e.sampleID).distinct().count();
			long mutCnt = get().getDataStream(e -> e.studyName.equals(study)).filter(SELECT_MUT::select).distinct().count();

			double average = mutCnt / (double) sampleCnt;

			System.out.println(study + "\t" + average + "\t" + sampleCnt);
		}
	}

	private static void printSampleOverlaps()
	{
		List<String> studies = get().getDataStream(e -> e.primaryPhenotype.equals("autism")).map(e -> e.studyName)
			.distinct().sorted().collect(Collectors.toList());

		for (String study1 : studies)
		{
			for (String study2 : studies)
			{
				if (study1.compareTo(study2) < 0)
				{
					Set<String> overlap = checkSampleOverlap(study1, study2);
					if (!overlap.isEmpty())
					{
						System.out.println(study1 + ", " + study2 + " = " + overlap);
					}
				}
			}
		}
	}

	private static Set<String> checkSampleOverlap(String study1, String study2)
	{
		Set<String> samples1 = get().getDataStream(e -> e.studyName.equals(study1)).map(e -> e.sampleID).collect(Collectors.toSet());
		Set<String> samples2 = get().getDataStream(e -> e.studyName.equals(study2)).map(e -> e.sampleID).collect(Collectors.toSet());

		return CollectionUtil.getIntersection(samples1, samples2);
	}

	private static void printMutationOverlapInCommonSamples(String study1, String study2)
	{
		Set<String> samples1 = get().getDataStream(e -> e.studyName.equals(study1)).map(e -> e.sampleID).collect(Collectors.toSet());
		Set<String> samples2 = get().getDataStream(e -> e.studyName.equals(study2)).map(e -> e.sampleID).collect(Collectors.toSet());

		Set<String> common = CollectionUtil.getIntersection(samples1, samples2);

		Set<String> muts1 = new HashSet<>();
		Set<String> muts2 = new HashSet<>();

		get().getDataStream(e -> (e.studyName.equals(study1) || e.studyName.equals(study2)) &&
			common.contains(e.sampleID) && SELECT_MUT.select(e)).forEach(e ->
		{
			String key = e.sampleID + " " + e.gene;
			Set<String> set = e.studyName.equals(study1) ? muts1 : muts2;
			set.add(key);
		});

		CollectionUtil.printNameMapping(study1, study2);
		CollectionUtil.printVennCounts(muts1, muts2);
	}

	private static void printRowsPerSamplePerStudy()
	{
		for (String study : GOOD_AUTISM_STUDIES)
		{
			Map<String, Long> sampleCounts = get().getDataStream(e -> e.studyName.equals(study) &&
				e.primaryPhenotype.equals("autism") &&
				(!e.studyName.equals("Yuen2017") || get().yuen2017T3SampleNames.contains(e.sampleID)))
				.collect(Collectors.groupingBy(e -> e.sampleID, Collectors.counting()));

			Histogram h = new Histogram(10);
			h.setBorderAtZero(true);
//			h.setUseLowerBorderForPrinting(true);
			sampleCounts.values().forEach(h::count);
			System.out.println("\nstudy = " + study);
			h.print();
		}
	}

	private static void printGenesPerSamplePerStudy()
	{
		DataFilter[] filters = new DataFilter[]{DataFilterEnum.YUEN_HAS_GENE, DataFilterEnum.AN_HAS_GENE};
		for (DataFilter filter : filters)
		{
			Map<String, List<Entry>> map = get().getDataStream(filter).collect(Collectors.groupingBy(e -> e.sampleID));

			Map<String, Long> counts = map.keySet().stream().collect(Collectors.toMap(Function.identity(),
				id -> map.get(id).stream().map(e -> e.gene).distinct().count()));

			Histogram h = new Histogram(5);
			h.setBorderAtZero(true);
			h.setUseLowerBorderForPrinting(true);
			counts.values().forEach(h::count);
			System.out.println("\nfilter = " + filter);
			h.print();
		}
	}

	private static void printMutationClassCounts()
	{
		DataFilter filter = e -> SELECT_MUT.select(e) && SELECT_AUTISM_STUDY.select(e) && e.gene.equals("NR3C1");
		Map<String, Long> cnt = get().getDataStream(filter).collect(Collectors.groupingBy(e -> e.functionClass, Collectors.counting()));
		cnt.keySet().stream().sorted(Comparator.comparing(cnt::get).reversed()).forEach(c -> System.out.println(cnt.get(c) + "\t" + c));

		long count = get().getDataStream(filter).count();
		System.out.println("\nTotal count = " + count);
	}

	private static void recordMutationNonUniformity() throws IOException
	{
		Map<String, Map<String, Integer>> map = new HashMap<>();
		get().getDataStream(e -> SELECT_AUTISM_STUDY.select(e) && SELECT_MUT.select(e)).forEach(e ->
		{
			if (!map.containsKey(e.gene)) map.put(e.gene, new HashMap<>());
			map.get(e.gene).put(e.sampleID, e.position);
		});

		Map<String, List<Integer>> posMap = map.keySet().stream().filter(gene -> map.get(gene).size() > 1)
			.filter(gene -> gene.equals("NR3C1"))
			.collect(Collectors.toMap(Function.identity(), gene -> new ArrayList<>(map.get(gene).values())));

		posMap.values().forEach(Collections::sort);

		Map<String, Double> pvals = new HashMap<>();
		Progress prg = new Progress(posMap.size(), "Calculating non-uniformity");
		posMap.keySet().forEach(gene ->
		{
			int[] loc = HG37.get().getLocation(gene);
			if (loc != null)
			{
				List<Integer> pos = posMap.get(gene);

				pos = pos.stream().filter(p -> p >= loc[0] && p <= loc[1]).collect(Collectors.toList());

//				UniformityChecker.plot(CollectionUtil.convertIntegerListToDouble(pos), min, max);

				if (pos.size() > 1)
				{
					double p = KolmogorovSmirnov.testUniformityOfLocations(pos, loc[0], loc[1]);
					pvals.put(gene, p);
				}
			}
			prg.tick();
		});

		BufferedWriter writer = Files.newBufferedWriter(Paths.get("/home/ozgun/Documents/Grants/Mental/non-uniformityXX.txt"));
		writer.write("Gene\tMutSize\tpval");
		pvals.keySet().stream().sorted(Comparator.comparing(pvals::get)).forEach(gene -> FileUtil.lnwrite(
			gene + "\t" + posMap.get(gene).size() + "\t" + pvals.get(gene), writer));
		writer.close();
	}

	private static void printMutationUniformityAroundAGene()
	{
		String gene = "NR3C1";

		String chr = HG37.get().getChromosome(gene);
		System.out.println("chr = " + chr);

		int[] loc = HG37.get().getLocation(gene);

		int diff = loc[1] - loc[0];

		int left = loc[0] - diff;
		int right = loc[1] + diff;

		Map<String, Integer> map = new HashMap<>();
		get().getDataStream(SELECT_AUTISM_STUDY).forEach(e ->
		{
			if (e.chr.equals(chr) && e.position > left && e.position < right)
			{
				map.put(e.sampleID, e.position);
			}
			else if (e.gene.equals(gene) && e.position >= loc[0] && e.position <= loc[1])
			{
				System.out.println("e.chr = " + e.chr);
			}
		});

		List<Integer> posList = new ArrayList<>(map.values());
		System.out.println("posList.size() = " + posList.size());

		posList.sort(Integer::compareTo);

		UniformityChecker.plot(CollectionUtil.convertIntegerListToDouble(posList), left, right);
	}

	private static void plotHotspotsForAGene()
	{
		String gene = "NR3C1";

		Map<String, Integer> map = new HashMap<>();
		get().getDataStream(SELECT_AUTISM_STUDY).filter(e -> e.gene.equals(gene)).forEach(e ->
			map.put(e.sampleID, e.position));

		List<Integer> posList = new ArrayList<>(map.values());
		System.out.println("posList.size() = " + posList.size());

		posList.sort(Integer::compareTo);

		double[] v = ArrayUtil.convertIntListToBasicDoubleArray(posList);

		List<Integer> rList = randomlyRedistribute(posList);
		double[] r = ArrayUtil.convertIntListToBasicDoubleArray(rList);

		KernelDensityPlot.plot(gene + " positions", v);
		KernelDensityPlot.plot(gene + " random positions", r);
	}

	private static List<Integer> randomlyRedistribute(List<Integer> orig)
	{
		int min = orig.stream().min(Integer::compareTo).get();
		int max = orig.stream().max(Integer::compareTo).get();

		int dif = max - min;

		Random r = new Random();
		List<Integer> list = new ArrayList<>();

		for (int i = 0; i < orig.size(); i++)
		{
			list.add(r.nextInt(dif + 1) + min);
		}

		return list;
	}

	private static void printStudiesContent()
	{
		DataFilter genOrEx = e -> e.sequenceType.equals("genome") ||e.sequenceType.equals("exome");
			Set<String> focus = new HashSet<>(Arrays.asList("autism", "bipolar_type1", "bipolar_type2", "schizophrenia"));
		Set<String> studies = get().getDataStream(e -> focus.contains(e.primaryPhenotype) && genOrEx.select(e)).map(e -> e.studyName).collect(Collectors.toSet());
		Map<String, Map<String, Set<String>>> phenCounts = new HashMap<>();
		Map<String, String> pubmedMap = new HashMap<>();
		Map<String, String> seqTypeMap = new HashMap<>();
		get().getDataStream(e -> studies.contains(e.studyName) && genOrEx.select(e)).forEach(e ->
		{
			if (!phenCounts.containsKey(e.studyName)) phenCounts.put(e.studyName, new HashMap<>());
			if (!phenCounts.get(e.studyName).containsKey(e.primaryPhenotype)) phenCounts.get(e.studyName).put(e.primaryPhenotype, new HashSet<>());
			phenCounts.get(e.studyName).get(e.primaryPhenotype).add(e.sampleID);
			pubmedMap.put(e.studyName, e.pubmedID);
			seqTypeMap.put(e.studyName, e.sequenceType);
		});


		phenCounts.forEach((study, map) ->
		{
			System.out.println("\n" + study + ", " + pubmedMap.get(study) + ", " + seqTypeMap.get(study));
			for (String phen : map.keySet())
			{
				System.out.println(phen + " = " + map.get(phen).size());
			}
		});
	}

	private static void temp()
	{
	}
}
