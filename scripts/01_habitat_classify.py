#!/usr/bin/env python3
"""Phase 2a: classify GTDB genomes with coordinates into habitat classes from
NCBI isolation_source text, parse coordinates, and flag depth/sediment.

Input : GTDB r232 bac120 metadata (filtered to CheckM2 >=90% / <5% with lat_lon)
Output: data_small/habitat_classified.tsv.gz

Rules are ordered by priority: host/clinical and engineered habitats are
assigned first so that e.g. "sea lion feces" is host-associated, not marine.
"""
import re, sys
import pandas as pd, numpy as np

SRC = sys.argv[1] if len(sys.argv) > 1 else '/home/user/data/gtdb/hq_latlon.tsv.gz'
OUT = sys.argv[2] if len(sys.argv) > 2 else 'data_small/habitat_classified.tsv.gz'

R = lambda p: re.compile(p, re.I)
RULES = [
    ('host_human_clinical', R(r'feces|faeces|fecal|faecal|stool|blood|sputum|urine|patient|clinical|hospital|wound|swab|pus\b|abscess|cerebrospinal|csf\b|biopsy|tissue|saliva|oral|dental|plaque|tongue|vagin|skin|nasal|nasophar|oropharyn|throat|bronch|lung|gut\b|intestin|colon|rectal|cecum|caecum|rumen|abomasum|stomach|gastric|bile|milk|udder|mastitis|infant|human|homo sapiens|bovine|cattle|chicken|poultry|pig\b|swine|porcine|mouse|mice|\brat\b|insect|termite|\bant\b|bee\b|honeybee|mosquito|tick\b|nematode|worm|sponge|coral|ascidian|tunicate|bryozoan|mussel|oyster|clam|shrimp|crab|fish|squid|whale|seal\b|sea lion|lion|zoo|animal|host|symbiont|endosymbiont|trophosome|gill|lichen|nodule|root\b|roots|leaf|leaves|phyllosphere|rhizoplane|endophyt|plant|seed|flower|fruit|stem|algae|alga\b|seaweed|kelp|diatom|cyanobacterial bloom|microcystis colony|phycosphere|carcass|egg|larva|cloaca|manure|dung|meconium|ear\b|eye\b|gland|secretion|aspirate|lavage|ixodes|tracheal|respirator|sinus|nose\b|rectum|urethra|ocular|caecal|cystic fibrosis|catheter|implant|joint|bone|heart|liver|kidney|spleen|lymph')),
    ('engineered', R(r'sludge|wastewater|waste water|sewage|effluent|influent|bioreactor|reactor|digester|anaerobic digest|biogas|compost|landfill|food|cheese|kimchi|fermented|ferment|brewery|beer|wine|vinegar|sausage|meat|pork|beef|dairy|yogurt|sauerkraut|pickle|silage|oil field|oilfield|oil reservoir|produced water|fracking|hydraulic fractur|mine tailing|tailings|e-waste|electronic waste|dumpsite|contaminated|polluted|pollut|industrial|factory|plant\b.*treatment|treatment plant|water treatment|drinking water|tap water|pipe|cooling tower|aquarium|aquaculture|biofloc|rearing|artificial sea|supermarket|farm worker|pond culture|culture collection|laboratory|lab\b|enrichment|bioremed|biofilter|biocathode|electrode|anode|cathode|fuel cell|spacecraft|clean ?room|international space station|hospital sink|sink|shower|swimming pool|subsea tunnel|tunnel|concrete|paint|petroleum|crude oil|diesel|kerosene|acid mine drainage|amd\b|mining|mine water|rees-amd')),
    ('hydrothermal', R(r'hydrothermal|hot spring|hotspring|geyser|thermal spring|vent\b|vents\b|chimney|fumarole|solfatar|geothermal|mud volcano')),
    ('subsurface', R(r'groundwater|ground water|aquifer|borehole|subsurface|deep biosphere|aspo|äspö|fracture fluid|mine\b|cave|karst|bedrock|well water|deep well|crustal|basalt|oceanic crust|coalbed|coal bed|shale|serpentin')),
    ('hypersaline', R(r'hypersaline|saltern|salt pan|salt lake|salt crust|brine|halite|solar salt|soda lake|salar|sabkha|evaporit|salt marsh|saline lake|saline soil|saline')),
    ('marine', R(r'envo:00002150|envo_00002150|seawater|sea water|ocean|marine|\bsea\b|sea\s|coastal|coast\b|estuar|gulf|lagoon|bay\b|fjord|mangrove|intertidal|tidal|seafloor|sea floor|seabed|abyss|hadal|trench|bathy|pelag|continental shelf|shelf|beach|salt water|saltwater|brackish|caspian|baltic|black sea|red sea|mediterranean|pacific|atlantic|arctic ocean|north sea|china sea|sargasso|offshore|cold seep|methane seep|seep\b|whale fall|sea ice|tara')),
    ('freshwater', R(r'freshwater|fresh water|lake|river|stream|creek|pond|reservoir|wetland|marsh|swamp|bog\b|fen\b|peatland|spring water|spring\b|glacier|glacial|meltwater|snow|ice\b|cryoconite|rain|puddle|waterfall|canal|ditch|paddy|rice field|inland water|limnetic|lentic|lotic')),
    ('soil', R(r'soil|permafrost|desert|dune|sand\b|rhizosphere|grassland|forest|farmland|field|agricultur|cropland|crop\b|tundra|steppe|savanna|meadow|pasture|orchard|vineyard|litter|humus|loess|clay|rock\b|rocks\b|stone|biocrust|crust\b|terrestrial|ground\b|dust|air\b|aerosol|atmospher|mud\b|dirt|woodland|wasteland|farm land|highland')),
]
SEDIMENT = R(r'sediment|mud\b|silt|ooze|sand\b|seafloor|sea floor|seabed|benthic|core\b|subseafloor')
DEEP = R(r'(?:deep|abyss|hadal|trench|bathy|mesopelag|abyssopelag|aphotic|bottom water|\b\d{3,5}\s?m\b|\b\d{3,5}\s?meters|seep|cold seep|whale fall|seamount)')
WATER = R(r'water|seawater|column|plankton|pelag|surface|fraction|filter|bloom|ice')


def parse_latlon(s):
    m = re.match(r'^\s*(-?\d+(?:\.\d+)?)\s*([NS])\s+(-?\d+(?:\.\d+)?)\s*([EW])', str(s))
    if not m:
        return np.nan, np.nan
    lat = float(m.group(1)) * (1 if m.group(2) == 'N' else -1)
    lon = float(m.group(3)) * (1 if m.group(4) == 'E' else -1)
    if abs(lat) > 90 or abs(lon) > 180:
        return np.nan, np.nan
    return lat, lon


def classify(src):
    s = str(src).lower().strip()
    if s in ('none', 'nan', '', 'missing', 'not applicable', 'not collected', 'unknown', 'na', 'n/a', 'environmental', 'environment', 'not provided', 'other'):
        return 'unknown'
    for name, rx in RULES:
        if rx.search(s):
            return name
    if 'water' in s:
        return 'water_unspecified'
    return 'other'


def main():
    m = pd.read_csv(SRC, sep='\t', low_memory=False)
    ll = m.ncbi_lat_lon.map(parse_latlon)
    m['lat'] = [x[0] for x in ll]
    m['lon'] = [x[1] for x in ll]
    m = m.dropna(subset=['lat', 'lon'])
    # drop the (0,0) placeholder coordinates
    m = m[~((m.lat == 0) & (m.lon == 0))]
    src = m.ncbi_isolation_source.astype(str)
    m['habitat'] = src.map(classify)
    m['is_sediment'] = src.str.contains(SEDIMENT)
    m['is_deep'] = src.str.contains(DEEP)
    m['is_water'] = src.str.contains(WATER)
    m['is_mag'] = m.ncbi_genome_category.astype(str).str.contains('metagenome|single cell|environmental')
    m['species'] = m.gtdb_taxonomy.str.split(';').str[-1].str.replace('s__', '', regex=False)
    m['phylum'] = m.gtdb_taxonomy.str.split(';').str[1].str.replace('p__', '', regex=False)
    m.to_csv(OUT, sep='\t', index=False, compression='gzip')
    print(m.habitat.value_counts().to_string())
    print('species per habitat:')
    print(m.groupby('habitat').species.nunique().to_string())
    print('marine water non-deep non-sediment species:', m[(m.habitat == 'marine') & ~m.is_deep & ~m.is_sediment].species.nunique())


if __name__ == '__main__':
    main()
