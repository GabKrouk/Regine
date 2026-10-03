// Transcription-factor family of a gene, inferred from its annotation (symbol, Araport11 short description,
// UniProt name; then the curator summary). First matching rule wins, so the order matters.
// The result is written once into data/tf_annotation.tsv, where it can be checked and corrected by hand.
export const FAMILIES = [
  ["HSF", /heat shock (transcription )?factor|\bHSF[ABC]?\d/i],
  ["NAC", /\bNAC\b|NAC domain|\bANAC\d|\bATAF\d|\bCUC\d|no apical meristem/i],
  ["WRKY", /WRKY/i],
  ["bZIP", /bZIP|basic.leucine.zipper|basic region.leucine|\bABF\d|\bHY5\b|\bTGA\d|\bGBF\d|\bVIP1\b|\bAREB\d|leucine zipper transcription factor/i],
  ["bHLH", /bHLH|basic helix.?loop.?helix/i],
  ["AP2/ERF", /AP2\/ERF|\bAP2\b|AP2 domain|\bERF\d*\b|\bDREB|ethylene.responsive element.binding|Integrase-type DNA-binding|ethylene response factor|\bRAP2|cytokinin response factor/i],
  ["B3", /\bB3\b|AP2\/B3|auxin response factor|\bARF\d|\bRAV\d|ABI3|VP1|\bREM\d|reproductive meristem/i],
  ["MADS", /MADS|\bAGL\d|AGAMOUS|SEPALLATA/i],
  ["TCP", /\bTCP\d*\b|TEOSINTE BRANCHED/i],
  ["SBP", /\bSBP\b|SQUAMOSA promoter.binding|\bSPL\d/i],
  ["LBD", /LOB domain|\bLBD\d|lateral organ boundar/i],
  ["GRAS", /\bGRAS\b/],
  ["Trihelix", /trihelix|\bGT-?\d\b|GT-2/i],
  ["E2F/DP", /\bE2F|\bDPB?\b|E2F.DP/i],
  ["ZF-HD", /ZF-HD|zinc finger homeodomain|homeobox protein 3[0-4]\b/i],
  ["MYB", /myb domain protein|\bMYB\d|R2R3|MYB3R/i],
  ["G2-like (GARP)", /\bPHL\d|\bHHO\d|\bMYR\d|\bKAN\d|KANADI|\bGLK\d|\bHRS1\b|\bNIGT|\bPHR1\b|GARP|G2-like|\bEFM\b/i],
  ["ARR-B", /\bARR\d+|response regulator/i],
  ["MYB-related", /MYB-related|myb-like|\bRVE\d|\bLHY\b|\bCCA1\b|\bLCL\d|telomer|\bTRF|\bTRB\d|\bTRP\d|homeodomain-like/i],
  ["HD (homeobox)", /homeobox|homeodomain(?!-like)|HD-ZIP|\bKNAT|\bBEL1|\bWOX\d|\bATHB/i],
  ["Dof", /\bDof\b|DOF\d|\bCDF\d/],
  ["GATA", /\bGATA\b|GATA\d|GATA type/i],
  ["YABBY", /YABBY/i],
  ["CO-like (BBX)", /CONSTANS|\bBBX\d|B-box/i],
  ["C2H2", /C2H2|\bZAT\d|\bIDD\d|indeterminate.?\(?ID\)?.domain/i],
  ["GRF", /growth.regulating factor|\bGRF\d/i],
  ["CPP", /\bCPP\b|tesmin|TSO1/i],
  ["EIL", /\bEIN3\b|\bEIL\d|ethylene insensitive 3/i],
  ["BES1", /\bBES1\b|\bBZR\d|\bBEH\d|\bBMY2\b/i],
  ["NF-Y", /NF-Y|nuclear factor Y|CCAAT/i],
  ["CAMTA", /CAMTA|calmodulin.binding transcription|calmodulin binding;transcription|calmodulin binding protein/i],
  ["RWP-RK", /RWP-RK|NIN.like|\bNLP\d/i],
  ["GeBP", /GeBP|GLABROUS1 enhancer/i],
  ["SRS", /\bSRS\d|SHI.RELATED|\bSTY\d/i],
  ["S1Fa-like", /S1FA/i],
  ["BBR-BPC", /\bBPC\d|basic pentacysteine|BBR\/BPC/i],
  ["FAR1", /\bFAR1\b|\bFRS\d|\bFHY3\b/i],
  ["HMG", /\bHMG|high mobility group/i],
  ["Whirly", /whirly/i],
  ["ARID", /\bARID\b/i],
  ["Storekeeper", /storekeeper|\bSTKL/i],
  ["PAH (SIN3-like)", /paired amphipathic helix|\bSNL\d/i],
  ["mTERF", /transcription termination factor/i],
  ["Zinc finger (other)", /zinc finger|zinc knuckle/i],
];

export function familyOf(agi, symbol, shortDesc, texts = []) {
  const prot = (texts[3] || "").split(" — ")[0];
  const first = [symbol, shortDesc, texts[0], prot].join(" | ");
  for (const [f, re] of FAMILIES) if (re.test(first)) return f;
  const rest = [texts[1], texts[2], texts[3]].join(" | ");
  for (const [f, re] of FAMILIES) if (re.test(rest)) return f;
  return "";
}
