#!/usr/bin/env python3

"""
Script to convert MochiView annotation file to GFF
Copyright 2018 Thaddeus D. Seher.

=== GFF3 file format ===

 Column 1: "seqid"
  The ID of the landmark used to establish the coordinate system for the current
  feature. IDs may contain any characters, but must escape any characters not
  in the set [a-zA-Z0-9.:^*$@!+_?-|]. In particular, IDs may not contain
  unescaped whitespace and must not begin with an unescaped ">".
 Column 2: "source"
  The source is a free text qualifier intended to describe the algorithm or
  operating procedure that generated this feature. Typically this is the name
  of a piece of software, such as "Genescan" or a database name, such as
  "Genbank." In effect, the source is used to extend the feature ontology by
  adding a qualifier to the type creating a new composite type that is a
  subclass of the type in the type column.
 Column 3: "type"
  The type of the feature (previously called the "method"). This is constrained
  to be either a term from the Sequence Ontology or an SO accession number. The
  latter alternative is distinguished using the syntax SO:000000. In either
  case, it must be sequence_feature (SO:0000110) or an is_a child of it.
 Columns 4 & 5: "start" and "end"
  The start and end coordinates of the feature are given in positive 1-based
  
  integer coordinates, relative to the landmark given in column one. Start is
  always less than or equal to end. For features that cross the origin of a
  circular feature (e.g. most bacterial genomes, plasmids, and some viral
  genomes), the requirement for start to be less than or equal to end is
  satisfied by making end = the position of the end + the length of the landmark
  feature.

  For zero-length features, such as insertion sites, start equals end and the
  implied site is to the right of the indicated base in the direction of the
  landmark.
 Column 6: "score"
  The score of the feature, a floating point number. As in earlier versions of
  the format, the semantics of the score are ill-defined. It is strongly
  recommended that E-values be used for sequence similarity features, and that
  P-values be used for ab initio gene prediction features.
 Column 7: "strand"
  The strand of the feature. + for positive strand (relative to the landmark),
  - for minus strand, and . for features that are not stranded. In addition,
  ? can be used for features whose strandedness is relevant, but unknown.
 Column 8: "phase"
  For features of type "CDS", the phase indicates where the feature begins with
  reference to the reading frame. The phase is one of the integers 0, 1, or 2,
  indicating the number of bases that should be removed from the beginning of
  this feature to reach the first base of the next codon. In other words, a
  phase of "0" indicates that the next codon begins at the first base of the
  region described by the current line, a phase of "1" indicates that the next
  codon begins at the second base of this region, and a phase of "2" indicates
  that the codon begins at the third base of this region. This is NOT to be
  confused with the frame, which is simply start modulo 3.

  For forward strand features, phase is counted from the start field. For
  reverse strand features, phase is counted from the end field.

  The phase is REQUIRED for all CDS features.

 Column 9: "attributes"
  A list of feature attributes in the format tag=value. Multiple tag=value pairs
  are separated by semicolons. URL escaping rules are used for tags or values
  containing the following characters: ",=;". Spaces are allowed in this field,
  but tabs must be replaced with the %09 URL escape. Attribute values do not
  need to be and should not be quoted. The quotes should be included as part of
  the value by parsers and not stripped.

  These tags have predefined meanings:

  ID
   Indicates the ID of the feature. The ID attribute is required for features
   that have children (e.g. gene and mRNAs), or for those that span multiple
   lines, but are optional for other features. IDs for each feature must be
   unique within the scope of the GFF file. In the case of discontinuous
   features (i.e. a single feature that exists over multiple genomic locations)
   the same ID may appear on multiple lines. All lines that share an ID must
   collectively represent a single feature.
  Name
   Display name for the feature. This is the name to be displayed to the user.
   Unlike IDs, there is no requirement that the Name be unique within the file.
  Alias
   A secondary name for the feature. It is suggested that this tag be used
   whenever a secondary identifier for the feature is needed, such as locus
   names and accession numbers. Unlike ID, there is no requirement that Alias
   be unique within the file.
  Parent
   Indicates the parent of the feature. A parent ID can be used to group exons
   into transcripts, transcripts into genes, an so forth. A feature may have
   multiple parents. Parent can only be used to indicate a partof relationship.
  Target
   Indicates the target of a nucleotide-to-nucleotide or protein-to-nucleotide
   alignment. The format of the value is "target_id start end [strand]", where
   strand is optional and may be "+" or "-". If the target_id contains spaces,
   they must be escaped as hex escape %20.
  Gap
   The alignment of the feature to the target if the two are not collinear
   (e.g. contain gaps). The alignment format is taken from the CIGAR format
   described in the Exonerate documentation. See "THE GAP ATTRIBUTE" for a
   description of this format.
  Derives_from
   Used to disambiguate the relationship between one feature and another when
   the relationship is a temporal one rather than a purely structural "part of"
   one. This is needed for polycistronic genes. See "PATHOLOGICAL CASES" for
   further discussion.
  Note
   A free text note.
  Dbxref
   A database cross reference. See the section "Ontology Associations and Db
   Cross References" for details on the format.
  Ontology_term
   A cross reference to an ontology term. See the section "Ontology Associations
   and Db Cross References" for details.
  Is_circular
   A flag to indicate whether a feature is circular. See extended discussion
   below.




=== MochiView file format ===
File Format:
 SEQ_NAME (required; max 50)
  The entry must correspond to the name of a Sequence in the database.
 START (required) and END (required)
  These columns must contain integer values, neither of which exceeds the length of the indicated Sequence.
  The 'START' and 'END' coordinates are inclusive, and the 'START' coordinate should be less than or equal to
  the 'END' coordinate.
 STRAND (optional)
  If this column is omitted, all Locations are assumed to be on the plus strand. If included, the value should be
  '+' (no quotes) for plus strand Locations and '-' for minus strand Locations (in practice, MochiView assigns
  anything without a minus to the plus strand).
 ANNO_TAG (optional; max 15) and ANNO_DESC (optional; max 500)
  These columns supply additional Location Annotations for the imported Locations. The value supplied for
  'ANNO_TAG' is displayed as an Annotation Tag on Locations when the Location Set is displayed on a Track.
  The value supplied for 'ANNO_DESC' is an Annotation Description that is visible when a Location in a plot is
  double-clicked. If only an 'ANNO_DESC' column is supplied (or the 'ANNO_TAG' column has a blank entry),
  the 'ANNO_TAG' is created from the first 15 characters of the 'ANNO_DESC' entry. Note that Location
  Annotations do not have to be supplied for all rows.
 GROUP (optional)
  If included, this column assigns Locations to groups. Currently the only effect of group information in
  MochiView is that Shape Style plots (described later) will connect the Locations in a group with a thick black
  line. This column must contain a numerical integer entry (above zero) for each row. Any Locations sharing
  the same integer are considered to be in the same group.

File Format: GENE FIELDS
 FEATURE_NAME (required; max 50)
  This column must contain a unique descriptor of the gene.
 GENE_NAME (optional v1.26+; max 50)
  The gene name should be unique if not left blank.
 ALIASES (optional; max 500)
  Aliases are provided as a pipe delimited text string of alternative names (e.g. HEX7|Contig12.5231|YUP7).
  Aliases do not have to be unique (but the program currently makes no use of non-unique aliases).
 DESCRIPTION (optional; max 500)
  This column can be left blank, or can contain a description of the gene that will be displayed in tooltips.

File Format: ISOFORM FIELDS
 TXN_START and TXN_END (required)
  These mandatory columns delineate the (inclusive) boundary coordinates of the isoform transcript. It does
  not matter which is the lower of the two values (MochiView infers the direction from the gene 'STRAND'
  entry). The coordinates for 'TXN_START' and 'TXN_END' must fall within the boundary delineated by the
  gene 'START' and 'END' coordinates.
 EXON_COUNT (required)
  This entry indicates the number of exons listed in the 'EXON_STARTS' and 'EXON_ENDS' fields.
 EXON_STARTS and EXON_ENDS (required)
  These two columns contain pipe-delimited coordinates for the coordinates (inclusive) of all exons in the
  isoform. As with the 'TXN_START' and 'TXN_END' columns, it does not matter which column contains the
  lower of the two values. For example, two exons with coordinates of 1000-->1500 and 2250-->3000 would be
  represented by an entry of '1000|2250' (omit quotes) in one column and '1500|3000' in the other (and an
  entry of '2' in the 'EXON_COUNT' column). The exon coordinates provided in these columns must all fall
  within the boundary delineated by the 'TXN_START' and 'TXN_END' columns.
 CDS_START and CDS_END (optional)
  These columns are optional, but if you include one you must include the other. The coordinates provided in
  these columns delineate the portion of the isoform transcript (if any) that is a protein coding region. These
  coordinates must fall within the boundary delineated by the 'EXON_STARTS ' and 'EXON_ENDS' columns.
  MochiView displays non-coding exons with higher transparency than coding exons. (If the boundary of the
  coding region falls within an exon, MochiView will split the exon.)
 ISOFORM_NAME (optional; max 50)
  Each isoform can be given a name (which can be seen in the plot tooltip for the isoform).
 IS_PRIMARY (optional)
  MochiView orders the isoforms in a vertically stacked display (using the criteria described below). An entry
  of 'Y' (omit quotes) assigns an isoform as 'primary', ensuring the isoform is displayed first (at the top of the
  stack). Only one isoform for each gene can have an entry of 'Y'.


type hierarchy suffix
 gene              (None)
  mRNA             -T
   exon            -T-E1, -T-E2, -T-E3, ...
   CDS             -P
   five_prime_UTR  -5
   three_prime_UTR -3

type list
 gene, mRNA, CDS, exon, intron, polyA_sequence, polyA_site, five_prime_UTR, three_prime_UTR, region, chromosome
"""

import sys
import os

__author__ = "Thaddeus D. Seher (@tdseher)"
__program__ = os.path.basename(sys.argv[0])

class feature_parser(object):
    def __init__(self, line, source):
        # Define MochiView column ID indices
        SEQ_NAME = 0
        START = 1
        END = 2
        STRAND = 3
        FEATURE_NAME = 4
        TXN_START = 5
        TXN_END = 6
        EXON_COUNT = 7
        EXON_STARTS = 8
        EXON_ENDS = 9
        CDS_START = 10
        CDS_END = 11
        GENE_NAME = 12
        ALIASES = 13
        DESCRIPTION = 14
        
        # Process the line
        line = line.rstrip()
        sline = line.split("\t")
        sline = [
            sline[SEQ_NAME],
            int(sline[START]),
            int(sline[END]),
            sline[STRAND],
            sline[FEATURE_NAME],
            int(sline[TXN_START]),
            int(sline[TXN_END]),
            int(sline[EXON_COUNT]),
            list(map(int, sline[EXON_STARTS].split('|'))),
            list(map(int, sline[EXON_ENDS].split('|'))),
            int(sline[CDS_START]) if (sline[CDS_START] != '') else None,
            int(sline[CDS_END]) if (sline[CDS_END] != '') else None,
            sline[GENE_NAME],
            sline[ALIASES],
            sline[DESCRIPTION],
        ]
        
        # Initialize feature list
        self.features = []
        
        base_attributes = {}
        if (sline[GENE_NAME] != ''):
            base_attributes['Name'] = sline[GENE_NAME]
        
        # Make sure CDS_START and CDS_END are in ascending order
        if (sline[CDS_START] != None):
            s_cds_start, s_cds_end = sorted([sline[CDS_START], sline[CDS_END]])
        # We assume the EXON_STARTS and EXON_ENDS are properly sorted
        # We make sure they are in ascending order also
        if (sline[STRAND] == '+'):
            s_exon_starts = sline[EXON_STARTS]
            s_exon_ends = sline[EXON_ENDS]
        elif (sline[STRAND] == '-'):
            s_exon_starts = sline[EXON_ENDS][::-1]
            s_exon_ends = sline[EXON_STARTS][::-1]
        
        ##### gene #####
        seqid = sline[SEQ_NAME]
        #source = 'MochiViewConverter'
        type = 'gene'
        start = sline[START]
        end = sline[END]
        score = '.'
        strand = sline[STRAND]
        phase = '.'
        attributes = base_attributes.copy()
        attributes['ID'] = sline[FEATURE_NAME]
        attributes['Alias'] = ','.join(sline[ALIASES].split('|'))
        attributes['Note'] = self._escape_text(sline[DESCRIPTION])
        
        self.features.append((seqid, source, type, start, end, score, strand, phase, self._join_attributes(attributes)))
        
        ##### mRNA #####
        type = 'mRNA' if (sline[CDS_START] != None) else 'RNA'
        start = sline[TXN_START]
        end = sline[TXN_END]
        score = '.'
        strand = sline[STRAND]
        phase = '.'
        attributes = base_attributes.copy()
        attributes['ID'] = sline[FEATURE_NAME]+'-T'
        attributes['Parent'] = sline[FEATURE_NAME]
        attributes['Alias'] = ','.join(sline[ALIASES].split('|'))
        attributes['Note'] = self._escape_text(sline[DESCRIPTION])
        
        self.features.append((seqid, source, type, start, end, score, strand, phase, self._join_attributes(attributes)))
        
        ##### exon #####
        assert int(sline[EXON_COUNT]) == len(sline[EXON_STARTS]) == len(sline[EXON_ENDS]), "Exons improperly defined"
        
        for i in range(sline[EXON_COUNT]):
            type = 'exon'
            start = s_exon_starts[i]
            end = s_exon_ends[i]
            score = '.'
            strand = sline[STRAND]
            phase = '.'
            attributes = base_attributes.copy()
            if (strand == '+'):
                attributes['ID'] = sline[FEATURE_NAME]+'-T-E'+str(i+1)
            elif (strand == '-'):
                attributes['ID'] = sline[FEATURE_NAME]+'-T-E'+str(sline[EXON_COUNT]-i)
            attributes['Parent'] = sline[FEATURE_NAME]+'-T'
            
            self.features.append((seqid, source, type, start, end, score, strand, phase, self._join_attributes(attributes)))
        
        # Add one CDS for every exon if it is a protein-coding gene
        if (sline[CDS_START] != None):
            if (sline[EXON_COUNT] == 1):
                type = 'CDS'
                start = s_cds_start
                end = s_cds_end
                score = '.'
                strand = sline[STRAND]
                phase = 0
                attributes = base_attributes.copy()
                attributes['ID'] = sline[FEATURE_NAME]+'-P'
                attributes['Parent'] = sline[FEATURE_NAME]+'-T'
                
                self.features.append((seqid, source, type, start, end, score, strand, phase, self._join_attributes(attributes)))
                
            #elif (int(sline[EXON_COUNT]) == 2):
            else:
                # Truncate the first exon if needed
                # Truncate the last exon if needed
                # If there are any in the middle, the keep the entire exon as a CDS
                
                # type = 'CDS'
                # start = s_cds_start
                # end = s_exon_ends[0]
                # score = '.'
                # strand = sline[STRAND]
                # phase = 0
                # attributes = base_attributes.copy()
                # attributes['ID'] = sline[FEATURE_NAME]+'-P'
                # attributes['Parent'] = sline[FEATURE_NAME]+'-T'
                # 
                # self.features.append((seqid, source, type, start, end, score, strand, phase, self._join_attributes(attributes)))
                
                for i in range(len(sline[EXON_STARTS])):
                    if ((s_cds_start < s_exon_ends[i]) and (s_exon_starts[i] < s_cds_end)):
                        
                        type = 'CDS'
                        
                        if (s_exon_starts[i] < s_cds_start < s_exon_ends[i]):
                            start = s_cds_start
                        else:
                            start = s_exon_starts[i]
                        
                        if (start < s_cds_end < s_exon_ends[i]):
                            end = s_cds_end
                        else:
                            end = s_exon_ends[i]
                        
                        #end = s_exon_ends[i]
                        score = '.'
                        strand = sline[STRAND]
                        phase = 0
                        attributes = base_attributes.copy()
                        attributes['ID'] = sline[FEATURE_NAME]+'-P'
                        attributes['Parent'] = sline[FEATURE_NAME]+'-T'
                        
                        self.features.append((seqid, source, type, start, end, score, strand, phase, self._join_attributes(attributes)))
                
                # Truncate the last exon
                # if (s_exon_starts[-1] < s_cds_end):
                #     type = 'CDS'
                #     start = s_exon_starts[-1]
                #     end = s_cds_end
                #     score = '.'
                #     strand = sline[STRAND]
                #     phase = 0
                #     attributes = base_attributes.copy()
                #     attributes['ID'] = sline[FEATURE_NAME]+'-P'
                #     attributes['Parent'] = sline[FEATURE_NAME]+'-T'
                #     
                #     self.features.append((seqid, source, type, start, end, score, strand, phase, self._join_attributes(attributes)))
            
            # Add upstream 5'(+)/3'(-) UTR
            if (sline[TXN_START] < s_cds_start):
                strand = sline[STRAND]
                start = sline[TXN_START]
                end = s_cds_start-1
                attributes = base_attributes.copy()
                if (strand == '+'):
                    type = 'five_prime_UTR'
                    #end = s_cds_start-1
                    attributes['ID'] = sline[FEATURE_NAME]+'-5'
                elif (strand == '-'):
                    type = 'three_prime_UTR'
                    #end = s_cds_start+1
                    attributes['ID'] = sline[FEATURE_NAME]+'-3'
                score = '.'
                phase = '.'
                attributes['Parent'] = sline[FEATURE_NAME]+'-T'
                
                self.features.append((seqid, source, type, start, end, score, strand, phase, self._join_attributes(attributes)))
            
            # Add downstream 3'(+)/5'(-) UTR
            if (s_cds_end < sline[TXN_END]):
                strand = sline[STRAND]
                start = s_cds_end+1
                end = sline[TXN_END]
                attributes = base_attributes.copy()
                if (strand == '+'):
                    type = 'three_prime_UTR'
                    #start = s_cds_end+1
                    attributes['ID'] = sline[FEATURE_NAME]+'-3'
                elif (strand == '-'):
                    type = 'five_prime_UTR'
                    #start = s_cds_end-1
                    attributes['ID'] = sline[FEATURE_NAME]+'-5'
                score = '.'
                phase = '.'
                attributes['Parent'] = sline[FEATURE_NAME]+'-T'
                
                self.features.append((seqid, source, type, start, end, score, strand, phase, self._join_attributes(attributes)))
    
    def _escape_text(self, text):
        '''
        GFF3 files are nine-column, tab-delimited, plain text files. Literal use
        of tab, newline, carriage return, the percent (%) sign, and control
        characters must be encoded using RFC 3986 Percent-Encoding; no other
        characters may be encoded. Backslash and other ad-hoc escaping
        conventions that have been added to the GFF format are not allowed.
        The file contents may include any character in the set supported by the
        operating environment, although for portability with other systems, use
        of Latin-1 or Unicode are recommended.
        
        Note that unescaped spaces are allowed within fields, meaning that
        parsers must split on tabs, not spaces. Use of the "+" (plus) character
        to encode spaces is depracated from early versions of the spec and is
        no longer allowed.
        
        Undefined fields are replaced with the "." character, as described in
        the original GFF spec.
        '''
        
        table = {
            '"': '%22',
            '%': '%25',
            '&': '%26',
            "'": '%27',
            '(': '%28',
            ')': '%29',
            ',': '%2C',
            ';': '%3B',
            '=': '%3D',
            '[': '%5B',
            '\\': '%5C',
            ']': '%5D',
            chr(127): '%7F',
        }
        for i in range(33):
            table[chr(i)] = '%' + hex(i)[2:].zfill(2)
        #table[chr(0)] = '%20'
        
        return ''.join(table.get(c,c) for c in text)
    
    def _join_attributes(self, attributes):
        l = []
        for k in attributes:
            l.append(str(k) + '=' + str(attributes[k]))
        return ';'.join(l)
    
    def __str__(self):
        return('\n'.join(map(lambda x: '\t'.join(map(str, x)), self.features)))

def main():
    if (len(sys.argv) != 3):
        print("USAGE: python3 {__program__} MochiViewAnnotations.txt source > MochiViewAnnotations.gff".format(**globals()))
        sys.exit()
    
    line_count = 0
    with open(sys.argv[1], 'r') as flo:
        for line in flo:
            if (line_count > 0):
                f = feature_parser(line, sys.argv[2]) # 'mochiview2gff-Tuch-et-al-2010'
                print(f)
            line_count += 1
    
    # Force encoding
    # with open('temp.new', 'w', encoding="ascii") as out:
    #     with open(sys.argv[1], 'r') as flo:
    #         for line in flo:
    #             if (line_count > 0):
    #                 f = feature_parser(line, sys.argv[2]) # 'mochiview2gff-Tuch-et-al-2010'
    #                 print(str(f), file=out)
    #             line_count += 1

if (__name__ == "__main__"):
    main()
