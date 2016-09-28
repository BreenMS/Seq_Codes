#!/usr/bin/perl -w
use strict;

my($file,$gi,$line,%hash, $out,$symbol, $gi2, $chr, $start, $stop,
   $geneid, $biotype, $probe,
   $s1, $s2, $s3, $s4, $s5, $s6, $s7, $s8, $s9, $s10, $s11, $s12, $s13, $s14, $s15, $s16, $s17, $s18, $s19, $s20, $s21, $s22, $s23, $s24, $s25, $s26, $s27, $s28, $s29, $s30, $s31, $s32, $s33, $s34, $s35, $s36, $s37, $s38, $s39, $s40, $s41, $s42, $s43, $s44,
   $s1p, $s2p, $s3p, $s4p, $s5p, $s6p, $s7p, $s8p, $s9p, $s10p, $s11p, $s12p, $s13p, $s14p, $s15p, $s16p, $s17p, $s18p, $s19p, $s20p, $s21p, $s22p, $s23p, $s24p, $s25p, $s26p, $s27p, $s28p, $s29p, $s30p, $s31p, $s32p, $s33p, $s34p, $s35p, $s36p, $s37p, $s38p, $s39p, $s40p, $s41p, $s42p, $s43p, $s44p,
   $s1ps, $s2ps, $s3ps, $s4ps, $s5ps, $s6ps, $s7ps, $s8ps, $s9ps, $s10ps, $s11ps, $s12ps, $s13ps, $s14ps, $s15ps, $s16ps, $s17ps, $s18ps, $s19ps, $s20ps, $s21ps, $s22ps, $s23ps, $s24ps, $s25ps, $s26ps, $s27ps, $s28ps, $s29ps, $s30ps, $s31ps, $s32ps, $s33ps, $s34ps, $s35ps, $s36ps, $s37ps, $s38ps, $s39ps, $s40ps, $s41ps, $s42ps, $s43ps, $s44ps,
   $s1nc, $s2nc, $s3nc, $s4nc, $s5nc, $s6nc, $s7nc, $s8nc, $s9nc, $s10nc, $s11nc, $s12nc, $s13nc, $s14nc, $s15nc, $s16nc, $s17nc, $s18nc, $s19nc, $s20nc, $s21nc, $s22nc, $s23nc, $s24nc, $s25nc, $s26nc, $s27nc, $s28nc, $s29nc, $s30nc, $s31nc, $s32nc, $s33nc, $s34nc, $s35nc, $s36nc, $s37nc, $s38nc, $s39nc, $s40nc, $s41nc, $s42nc, $s43nc, $s44nc,
   $s1tcp, $s2tcp, $s3tcp, $s4tcp, $s5tcp, $s6tcp, $s7tcp, $s8tcp, $s9tcp, $s10tcp, $s11tcp, $s12tcp, $s13tcp, $s14tcp, $s15tcp, $s16tcp, $s17tcp, $s18tcp, $s19tcp, $s20tcp, $s21tcp, $s22tcp, $s23tcp, $s24tcp, $s25tcp, $s26tcp, $s27tcp, $s28tcp, $s29tcp, $s30tcp, $s31tcp, $s32tcp, $s33tcp, $s34tcp, $s35tcp, $s36tcp, $s37tcp, $s38tcp, $s39tcp, $s40tcp, $s41tcp, $s42tcp, $s43tcp, $s44tcp,
   $s1probe, $s2probe, $s3probe, $s4probe, $s5probe, $s6probe, $s7probe, $s8probe, $s9probe, $s10probe, $s11probe, $s12probe, $s13probe, $s14probe, $s15probe, $s16probe, $s17probe, $s18probe, $s19probe, $s20probe, $s21probe, $s22probe, $s23probe, $s24probe, $s25probe, $s26probe, $s27probe, $s28probe, $s29probe, $s30probe, $s31probe, $s32probe, $s33probe, $s34probe, $s35probe, $s36probe, $s37probe, $s38probe, $s39probe, $s40probe, $s41probe, $s42probe, $s43probe, $s44probe,
   $s1pprobe, $s2pprobe, $s3pprobe, $s4pprobe, $s5pprobe, $s6pprobe, $s7pprobe, $s8pprobe, $s9pprobe, $s10pprobe, $s11pprobe, $s12pprobe, $s13pprobe, $s14pprobe, $s15pprobe, $s16pprobe, $s17pprobe, $s18pprobe, $s19pprobe, $s20pprobe, $s21pprobe, $s22pprobe, $s23pprobe, $s24pprobe, $s25pprobe, $s26pprobe, $s27pprobe, $s28pprobe, $s29pprobe, $s30pprobe, $s31pprobe, $s32pprobe, $s33pprobe, $s34pprobe, $s35pprobe, $s36pprobe, $s37pprobe, $s38pprobe, $s39pprobe, $s40pprobe, $s41pprobe, $s42pprobe, $s43pprobe, $s44pprobe,
   $s1psprobe, $s2psprobe, $s3psprobe, $s4psprobe, $s5psprobe, $s6psprobe, $s7psprobe, $s8psprobe, $s9psprobe, $s10psprobe, $s11psprobe, $s12psprobe, $s13psprobe, $s14psprobe, $s15psprobe, $s16psprobe, $s17psprobe, $s18psprobe, $s19psprobe, $s20psprobe, $s21psprobe, $s22psprobe, $s23psprobe, $s24psprobe, $s25psprobe, $s26psprobe, $s27psprobe, $s28psprobe, $s29psprobe, $s30psprobe, $s31psprobe, $s32psprobe, $s33psprobe, $s34psprobe, $s35psprobe, $s36psprobe, $s37psprobe, $s38psprobe, $s39psprobe, $s40psprobe, $s41psprobe, $s42psprobe, $s43psprobe, $s44psprobe,
   $s1ncprobe, $s2ncprobe, $s3ncprobe, $s4ncprobe, $s5ncprobe, $s6ncprobe, $s7ncprobe, $s8ncprobe, $s9ncprobe, $s10ncprobe, $s11ncprobe, $s12ncprobe, $s13ncprobe, $s14ncprobe, $s15ncprobe, $s16ncprobe, $s17ncprobe, $s18ncprobe, $s19ncprobe, $s20ncprobe, $s21ncprobe, $s22ncprobe, $s23ncprobe, $s24ncprobe, $s25ncprobe, $s26ncprobe, $s27ncprobe, $s28ncprobe, $s29ncprobe, $s30ncprobe, $s31ncprobe, $s32ncprobe, $s33ncprobe, $s34ncprobe, $s35ncprobe, $s36ncprobe, $s37ncprobe, $s38ncprobe, $s39ncprobe, $s40ncprobe, $s41ncprobe, $s42ncprobe, $s43ncprobe, $s44ncprobe,
   $s1tcpprobe, $s2tcpprobe, $s3tcpprobe, $s4tcpprobe, $s5tcpprobe, $s6tcpprobe, $s7tcpprobe, $s8tcpprobe, $s9tcpprobe, $s10tcpprobe, $s11tcpprobe, $s12tcpprobe, $s13tcpprobe, $s14tcpprobe, $s15tcpprobe, $s16tcpprobe, $s17tcpprobe, $s18tcpprobe, $s19tcpprobe, $s20tcpprobe, $s21tcpprobe, $s22tcpprobe, $s23tcpprobe, $s24tcpprobe, $s25tcpprobe, $s26tcpprobe, $s27tcpprobe, $s28tcpprobe, $s29tcpprobe, $s30tcpprobe, $s31tcpprobe, $s32tcpprobe, $s33tcpprobe, $s34tcpprobe, $s35tcpprobe, $s36tcpprobe, $s37tcpprobe, $s38tcpprobe, $s39tcpprobe, $s40tcpprobe, $s41tcpprobe, $s42tcpprobe, $s43tcpprobe, $s44tcpprobe);

open(INFILE, $ARGV[0]) or die"File1 is Dead\n";
while(<INFILE>)
{
$line=$_;
chomp $line;
if ($line=~m/\S+/)
{
    ($geneid)=  ($line=~ m/(\S+)/);
    ($biotype)= ($line=~m/\S+\s+(\S+)/);
    ($probe)=   ($line=~m/\S+\s+\S+\s+(\S+)/);

    
   if ($biotype=~/^lncRNA$/ or $biotype=~/^lncRNA_RefSeq$/ or $biotype=~/^lncRNA_Cabili$/)
   {
    # split the line into two fields                                                                                                    
    my @columns = split(/\t/, $line);

    if ($columns[3] > 1) {$s1++;}
    if ($columns[4] > 1) { $s2++;}
    if ($columns[5] > 1) { $s3++;}
    if ($columns[6] > 1) { $s4++;}
    if ($columns[7] > 1) { $s5++;}
    if ($columns[8] > 1) { $s6++;}
    if ($columns[9] > 1) { $s7++;}
    if ($columns[10] > 1) { $s8++;}
    if ($columns[11] > 1) { $s9++;}
    if ($columns[12] > 1) { $s10++;}
    if ($columns[13] > 1) { $s11++;}
    if ($columns[14] > 1) { $s12++;}
    if ($columns[15] > 1) { $s13++;}
    if ($columns[16] > 1) { $s14++;}
    if ($columns[17] > 1) { $s15++;}
    if ($columns[18] > 1) { $s16++;}
    if ($columns[19] > 1) { $s17++;}
    if ($columns[20] > 1) { $s18++;}
    if ($columns[21] > 1) { $s19++;}
    if ($columns[22] > 1) { $s20++;}
    if ($columns[23] > 1) { $s21++;}
    if ($columns[24] > 1) { $s22++;}
    if ($columns[25] > 1) { $s23++;}
    if ($columns[26] > 1) { $s24++;}
    if ($columns[27] > 1) { $s25++;}
    if ($columns[28] > 1) { $s26++;}
    if ($columns[29] > 1) { $s27++;}
    if ($columns[30] > 1) { $s28++;}
    if ($columns[31] > 1) { $s29++;}
    if ($columns[32] > 1) { $s30++;}
    if ($columns[33] > 1) { $s31++;}
    if ($columns[34] > 1) { $s32++;}
    if ($columns[35] > 1) { $s33++;}
    if ($columns[36] > 1) { $s34++;}
    if ($columns[37] > 1) { $s35++;}
    if ($columns[38] > 1) { $s36++;}
    if ($columns[39] > 1) { $s37++;}
    if ($columns[40] > 1) { $s38++;}
    if ($columns[41] > 1) { $s39++;}
    if ($columns[42] > 1) { $s40++;}
        if ($columns[43] > 1) { $s41++;}
        if ($columns[44] > 1) { $s42++;}
        if ($columns[45] > 1) { $s43++;}
        if ($columns[46] > 1) { $s44++;}
}
       if ($biotype=~/^protein_coding/)
   {
    # split the line into two fields                                                                                                    
    my @columns = split(/\t/, $line);

    if ($columns[3] > 1) {$s1p++;}
    if ($columns[4] > 1) { $s2p++;}
    if ($columns[5] > 1) { $s3p++;}
    if ($columns[6] > 1) { $s4p++;}
    if ($columns[7] > 1) { $s5p++;}
    if ($columns[8] > 1) { $s6p++;}
    if ($columns[9] > 1) { $s7p++;}
    if ($columns[10] > 1) { $s8p++;}
    if ($columns[11] > 1) { $s9p++;}
    if ($columns[12] > 1) { $s10p++;}
    if ($columns[13] > 1) { $s11p++;}
    if ($columns[14] > 1) { $s12p++;}
    if ($columns[15] > 1) { $s13p++;}
    if ($columns[16] > 1) { $s14p++;}
    if ($columns[17] > 1) { $s15p++;}
    if ($columns[18] > 1) { $s16p++;}
    if ($columns[19] > 1) { $s17p++;}
    if ($columns[20] > 1) { $s18p++;}
    if ($columns[21] > 1) { $s19p++;}
    if ($columns[22] > 1) { $s20p++;}
    if ($columns[23] > 1) { $s21p++;}
    if ($columns[24] > 1) { $s22p++;}
    if ($columns[25] > 1) { $s23p++;}
    if ($columns[26] > 1) { $s24p++;}
    if ($columns[27] > 1) { $s25p++;}
    if ($columns[28] > 1) { $s26p++;}
    if ($columns[29] > 1) { $s27p++;}
    if ($columns[30] > 1) { $s28p++;}
    if ($columns[31] > 1) { $s29p++;}
    if ($columns[32] > 1) { $s30p++;}
    if ($columns[33] > 1) { $s31p++;}
    if ($columns[34] > 1) { $s32p++;}
    if ($columns[35] > 1) { $s33p++;}
    if ($columns[36] > 1) { $s34p++;}
    if ($columns[37] > 1) { $s35p++;}
    if ($columns[38] > 1) { $s36p++;}
    if ($columns[39] > 1) { $s37p++;}
    if ($columns[40] > 1) { $s38p++;}
    if ($columns[41] > 1) { $s39p++;}
    if ($columns[42] > 1) { $s40p++;}
        if ($columns[43] > 1) { $s41p++;}
        if ($columns[44] > 1) { $s42p++;}
        if ($columns[45] > 1) { $s43p++;}
        if ($columns[46] > 1) { $s44p++;}
}
    
    
    
           if ($biotype=~/^pseudogene/)
   {
    # split the line into two fields                                                                                                    
    my @columns = split(/\t/, $line);

    if ($columns[3] > 1) {$s1ps++;}
    if ($columns[4] > 1) { $s2ps++;}
    if ($columns[5] > 1) { $s3ps++;}
    if ($columns[6] > 1) { $s4ps++;}
    if ($columns[7] > 1) { $s5ps++;}
    if ($columns[8] > 1) { $s6ps++;}
    if ($columns[9] > 1) { $s7ps++;}
    if ($columns[10] > 1) { $s8ps++;}
    if ($columns[11] > 1) { $s9ps++;}
    if ($columns[12] > 1) { $s10ps++;}
    if ($columns[13] > 1) { $s11ps++;}
    if ($columns[14] > 1) { $s12ps++;}
    if ($columns[15] > 1) { $s13ps++;}
    if ($columns[16] > 1) { $s14ps++;}
    if ($columns[17] > 1) { $s15ps++;}
    if ($columns[18] > 1) { $s16ps++;}
    if ($columns[19] > 1) { $s17ps++;}
    if ($columns[20] > 1) { $s18ps++;}
    if ($columns[21] > 1) { $s19ps++;}
    if ($columns[22] > 1) { $s20ps++;}
    if ($columns[23] > 1) { $s21ps++;}
    if ($columns[24] > 1) { $s22ps++;}
    if ($columns[25] > 1) { $s23ps++;}
    if ($columns[26] > 1) { $s24ps++;}
    if ($columns[27] > 1) { $s25ps++;}
    if ($columns[28] > 1) { $s26ps++;}
    if ($columns[29] > 1) { $s27ps++;}
    if ($columns[30] > 1) { $s28ps++;}
    if ($columns[31] > 1) { $s29ps++;}
    if ($columns[32] > 1) { $s30ps++;}
    if ($columns[33] > 1) { $s31ps++;}
    if ($columns[34] > 1) { $s32ps++;}
    if ($columns[35] > 1) { $s33ps++;}
    if ($columns[36] > 1) { $s34ps++;}
    if ($columns[37] > 1) { $s35ps++;}
    if ($columns[38] > 1) { $s36ps++;}
    if ($columns[39] > 1) { $s37ps++;}
    if ($columns[40] > 1) { $s38ps++;}
    if ($columns[41] > 1) { $s39ps++;}
    if ($columns[42] > 1) { $s40ps++;}
        if ($columns[43] > 1) { $s41ps++;}
        if ($columns[44] > 1) { $s42ps++;}
        if ($columns[45] > 1) { $s43ps++;}
        if ($columns[46] > 1) { $s44ps++;}
}
   
              if ($biotype=~/^ncRNA/)
   {
    # split the line into two fields                                                                                                    
    my @columns = split(/\t/, $line);

    if ($columns[3] > 1) {$s1nc++;}
    if ($columns[4] > 1) { $s2nc++;}
    if ($columns[5] > 1) { $s3nc++;}
    if ($columns[6] > 1) { $s4nc++;}
    if ($columns[7] > 1) { $s5nc++;}
    if ($columns[8] > 1) { $s6nc++;}
    if ($columns[9] > 1) { $s7nc++;}
    if ($columns[10] > 1) { $s8nc++;}
    if ($columns[11] > 1) { $s9nc++;}
    if ($columns[12] > 1) { $s10nc++;}
    if ($columns[13] > 1) { $s11nc++;}
    if ($columns[14] > 1) { $s12nc++;}
    if ($columns[15] > 1) { $s13nc++;}
    if ($columns[16] > 1) { $s14nc++;}
    if ($columns[17] > 1) { $s15nc++;}
    if ($columns[18] > 1) { $s16nc++;}
    if ($columns[19] > 1) { $s17nc++;}
    if ($columns[20] > 1) { $s18nc++;}
    if ($columns[21] > 1) { $s19nc++;}
    if ($columns[22] > 1) { $s20nc++;}
    if ($columns[23] > 1) { $s21nc++;}
    if ($columns[24] > 1) { $s22nc++;}
    if ($columns[25] > 1) { $s23nc++;}
    if ($columns[26] > 1) { $s24nc++;}
    if ($columns[27] > 1) { $s25nc++;}
    if ($columns[28] > 1) { $s26nc++;}
    if ($columns[29] > 1) { $s27nc++;}
    if ($columns[30] > 1) { $s28nc++;}
    if ($columns[31] > 1) { $s29nc++;}
    if ($columns[32] > 1) { $s30nc++;}
    if ($columns[33] > 1) { $s31nc++;}
    if ($columns[34] > 1) { $s32nc++;}
    if ($columns[35] > 1) { $s33nc++;}
    if ($columns[36] > 1) { $s34nc++;}
    if ($columns[37] > 1) { $s35nc++;}
    if ($columns[38] > 1) { $s36nc++;}
    if ($columns[39] > 1) { $s37nc++;}
    if ($columns[40] > 1) { $s38nc++;}
    if ($columns[41] > 1) { $s39nc++;}
    if ($columns[42] > 1) { $s40nc++;}
        if ($columns[43] > 1) { $s41nc++;}
        if ($columns[44] > 1) { $s42nc++;}
        if ($columns[45] > 1) { $s43nc++;}
        if ($columns[46] > 1) { $s44nc++;}
}
   
   
                 if ($biotype=~/^TUCP/)
   {
    # split the line into two fields                                                                                                    
    my @columns = split(/\t/, $line);

    if ($columns[3] > 1) {$s1tcp++;}
    if ($columns[4] > 1) { $s2tcp++;}
    if ($columns[5] > 1) { $s3tcp++;}
    if ($columns[6] > 1) { $s4tcp++;}
    if ($columns[7] > 1) { $s5tcp++;}
    if ($columns[8] > 1) { $s6tcp++;}
    if ($columns[9] > 1) { $s7tcp++;}
    if ($columns[10] > 1) { $s8tcp++;}
    if ($columns[11] > 1) { $s9tcp++;}
    if ($columns[12] > 1) { $s10tcp++;}
    if ($columns[13] > 1) { $s11tcp++;}
    if ($columns[14] > 1) { $s12tcp++;}
    if ($columns[15] > 1) { $s13tcp++;}
    if ($columns[16] > 1) { $s14tcp++;}
    if ($columns[17] > 1) { $s15tcp++;}
    if ($columns[18] > 1) { $s16tcp++;}
    if ($columns[19] > 1) { $s17tcp++;}
    if ($columns[20] > 1) { $s18tcp++;}
    if ($columns[21] > 1) { $s19tcp++;}
    if ($columns[22] > 1) { $s20tcp++;}
    if ($columns[23] > 1) { $s21tcp++;}
    if ($columns[24] > 1) { $s22tcp++;}
    if ($columns[25] > 1) { $s23tcp++;}
    if ($columns[26] > 1) { $s24tcp++;}
    if ($columns[27] > 1) { $s25tcp++;}
    if ($columns[28] > 1) { $s26tcp++;}
    if ($columns[29] > 1) { $s27tcp++;}
    if ($columns[30] > 1) { $s28tcp++;}
    if ($columns[31] > 1) { $s29tcp++;}
    if ($columns[32] > 1) { $s30tcp++;}
    if ($columns[33] > 1) { $s31tcp++;}
    if ($columns[34] > 1) { $s32tcp++;}
    if ($columns[35] > 1) { $s33tcp++;}
    if ($columns[36] > 1) { $s34tcp++;}
    if ($columns[37] > 1) { $s35tcp++;}
    if ($columns[38] > 1) { $s36tcp++;}
    if ($columns[39] > 1) { $s37tcp++;}
    if ($columns[40] > 1) { $s38tcp++;}
    if ($columns[41] > 1) { $s39tcp++;}
    if ($columns[42] > 1) { $s40tcp++;}
        if ($columns[43] > 1) { $s41tcp++;}
        if ($columns[44] > 1) { $s42tcp++;}
        if ($columns[45] > 1) { $s43tcp++;}
        if ($columns[46] > 1) { $s44tcp++;}
}
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
      if ($biotype=~/^lncRNA$/ or $biotype=~/^lncRNA_RefSeq$/ or $biotype=~/^lncRNA_Cabili$/ and $probe =~/^Y$/)
   {
    # split the line into two fields                                                                                                    
    my @columns = split(/\t/, $line);

    if ($columns[3] > 1) {$s1probe++;}
    if ($columns[4] > 1) { $s2probe++;}
    if ($columns[5] > 1) { $s3probe++;}
    if ($columns[6] > 1) { $s4probe++;}
    if ($columns[7] > 1) { $s5probe++;}
    if ($columns[8] > 1) { $s6probe++;}
    if ($columns[9] > 1) { $s7probe++;}
    if ($columns[10] > 1) { $s8probe++;}
    if ($columns[11] > 1) { $s9probe++;}
    if ($columns[12] > 1) { $s10probe++;}
    if ($columns[13] > 1) { $s11probe++;}
    if ($columns[14] > 1) { $s12probe++;}
    if ($columns[15] > 1) { $s13probe++;}
    if ($columns[16] > 1) { $s14probe++;}
    if ($columns[17] > 1) { $s15probe++;}
    if ($columns[18] > 1) { $s16probe++;}
    if ($columns[19] > 1) { $s17probe++;}
    if ($columns[20] > 1) { $s18probe++;}
    if ($columns[21] > 1) { $s19probe++;}
    if ($columns[22] > 1) { $s20probe++;}
    if ($columns[23] > 1) { $s21probe++;}
    if ($columns[24] > 1) { $s22probe++;}
    if ($columns[25] > 1) { $s23probe++;}
    if ($columns[26] > 1) { $s24probe++;}
    if ($columns[27] > 1) { $s25probe++;}
    if ($columns[28] > 1) { $s26probe++;}
    if ($columns[29] > 1) { $s27probe++;}
    if ($columns[30] > 1) { $s28probe++;}
    if ($columns[31] > 1) { $s29probe++;}
    if ($columns[32] > 1) { $s30probe++;}
    if ($columns[33] > 1) { $s31probe++;}
    if ($columns[34] > 1) { $s32probe++;}
    if ($columns[35] > 1) { $s33probe++;}
    if ($columns[36] > 1) { $s34probe++;}
    if ($columns[37] > 1) { $s35probe++;}
    if ($columns[38] > 1) { $s36probe++;}
    if ($columns[39] > 1) { $s37probe++;}
    if ($columns[40] > 1) { $s38probe++;}
    if ($columns[41] > 1) { $s39probe++;}
    if ($columns[42] > 1) { $s40probe++;}
        if ($columns[43] > 1) { $s41probe++;}
        if ($columns[44] > 1) { $s42probe++;}
        if ($columns[45] > 1) { $s43probe++;}
        if ($columns[46] > 1) { $s44probe++;}
}
       if ($biotype=~/^protein_coding/ and $probe =~/^Y$/)
   {
    # split the line into two fields                                                                                                    
    my @columns = split(/\t/, $line);

    if ($columns[3] > 1) {$s1pprobe++;}
    if ($columns[4] > 1) { $s2pprobe++;}
    if ($columns[5] > 1) { $s3pprobe++;}
    if ($columns[6] > 1) { $s4pprobe++;}
    if ($columns[7] > 1) { $s5pprobe++;}
    if ($columns[8] > 1) { $s6pprobe++;}
    if ($columns[9] > 1) { $s7pprobe++;}
    if ($columns[10] > 1) { $s8pprobe++;}
    if ($columns[11] > 1) { $s9pprobe++;}
    if ($columns[12] > 1) { $s10pprobe++;}
    if ($columns[13] > 1) { $s11pprobe++;}
    if ($columns[14] > 1) { $s12pprobe++;}
    if ($columns[15] > 1) { $s13pprobe++;}
    if ($columns[16] > 1) { $s14pprobe++;}
    if ($columns[17] > 1) { $s15pprobe++;}
    if ($columns[18] > 1) { $s16pprobe++;}
    if ($columns[19] > 1) { $s17pprobe++;}
    if ($columns[20] > 1) { $s18pprobe++;}
    if ($columns[21] > 1) { $s19pprobe++;}
    if ($columns[22] > 1) { $s20pprobe++;}
    if ($columns[23] > 1) { $s21pprobe++;}
    if ($columns[24] > 1) { $s22pprobe++;}
    if ($columns[25] > 1) { $s23pprobe++;}
    if ($columns[26] > 1) { $s24pprobe++;}
    if ($columns[27] > 1) { $s25pprobe++;}
    if ($columns[28] > 1) { $s26pprobe++;}
    if ($columns[29] > 1) { $s27pprobe++;}
    if ($columns[30] > 1) { $s28pprobe++;}
    if ($columns[31] > 1) { $s29pprobe++;}
    if ($columns[32] > 1) { $s30pprobe++;}
    if ($columns[33] > 1) { $s31pprobe++;}
    if ($columns[34] > 1) { $s32pprobe++;}
    if ($columns[35] > 1) { $s33pprobe++;}
    if ($columns[36] > 1) { $s34pprobe++;}
    if ($columns[37] > 1) { $s35pprobe++;}
    if ($columns[38] > 1) { $s36pprobe++;}
    if ($columns[39] > 1) { $s37pprobe++;}
    if ($columns[40] > 1) { $s38pprobe++;}
    if ($columns[41] > 1) { $s39pprobe++;}
    if ($columns[42] > 1) { $s40pprobe++;}
        if ($columns[43] > 1) { $s41pprobe++;}
        if ($columns[44] > 1) { $s42pprobe++;}
        if ($columns[45] > 1) { $s43pprobe++;}
        if ($columns[46] > 1) { $s44pprobe++;}
}
    
    
    
           if ($biotype=~/^pseudogene/ and $probe =~/^Y$/)
   {
    # split the line into two fields                                                                                                    
    my @columns = split(/\t/, $line);

    if ($columns[3] > 1) {$s1psprobe++;}
    if ($columns[4] > 1) { $s2psprobe++;}
    if ($columns[5] > 1) { $s3psprobe++;}
    if ($columns[6] > 1) { $s4psprobe++;}
    if ($columns[7] > 1) { $s5psprobe++;}
    if ($columns[8] > 1) { $s6psprobe++;}
    if ($columns[9] > 1) { $s7psprobe++;}
    if ($columns[10] > 1) { $s8psprobe++;}
    if ($columns[11] > 1) { $s9psprobe++;}
    if ($columns[12] > 1) { $s10psprobe++;}
    if ($columns[13] > 1) { $s11psprobe++;}
    if ($columns[14] > 1) { $s12psprobe++;}
    if ($columns[15] > 1) { $s13psprobe++;}
    if ($columns[16] > 1) { $s14psprobe++;}
    if ($columns[17] > 1) { $s15psprobe++;}
    if ($columns[18] > 1) { $s16psprobe++;}
    if ($columns[19] > 1) { $s17psprobe++;}
    if ($columns[20] > 1) { $s18psprobe++;}
    if ($columns[21] > 1) { $s19psprobe++;}
    if ($columns[22] > 1) { $s20psprobe++;}
    if ($columns[23] > 1) { $s21psprobe++;}
    if ($columns[24] > 1) { $s22psprobe++;}
    if ($columns[25] > 1) { $s23psprobe++;}
    if ($columns[26] > 1) { $s24psprobe++;}
    if ($columns[27] > 1) { $s25psprobe++;}
    if ($columns[28] > 1) { $s26psprobe++;}
    if ($columns[29] > 1) { $s27psprobe++;}
    if ($columns[30] > 1) { $s28psprobe++;}
    if ($columns[31] > 1) { $s29psprobe++;}
    if ($columns[32] > 1) { $s30psprobe++;}
    if ($columns[33] > 1) { $s31psprobe++;}
    if ($columns[34] > 1) { $s32psprobe++;}
    if ($columns[35] > 1) { $s33psprobe++;}
    if ($columns[36] > 1) { $s34psprobe++;}
    if ($columns[37] > 1) { $s35psprobe++;}
    if ($columns[38] > 1) { $s36psprobe++;}
    if ($columns[39] > 1) { $s37psprobe++;}
    if ($columns[40] > 1) { $s38psprobe++;}
    if ($columns[41] > 1) { $s39psprobe++;}
    if ($columns[42] > 1) { $s40psprobe++;}
        if ($columns[43] > 1) { $s41psprobe++;}
        if ($columns[44] > 1) { $s42psprobe++;}
        if ($columns[45] > 1) { $s43psprobe++;}
        if ($columns[46] > 1) { $s44psprobe++;}
}
   
              if ($biotype=~/^ncRNA/ and $probe =~/^Y$/)
   {
    # split the line into two fields                                                                                                    
    my @columns = split(/\t/, $line);

    if ($columns[3] > 1) {$s1ncprobe++;}
    if ($columns[4] > 1) { $s2ncprobe++;}
    if ($columns[5] > 1) { $s3ncprobe++;}
    if ($columns[6] > 1) { $s4ncprobe++;}
    if ($columns[7] > 1) { $s5ncprobe++;}
    if ($columns[8] > 1) { $s6ncprobe++;}
    if ($columns[9] > 1) { $s7ncprobe++;}
    if ($columns[10] > 1) { $s8ncprobe++;}
    if ($columns[11] > 1) { $s9ncprobe++;}
    if ($columns[12] > 1) { $s10ncprobe++;}
    if ($columns[13] > 1) { $s11ncprobe++;}
    if ($columns[14] > 1) { $s12ncprobe++;}
    if ($columns[15] > 1) { $s13ncprobe++;}
    if ($columns[16] > 1) { $s14ncprobe++;}
    if ($columns[17] > 1) { $s15ncprobe++;}
    if ($columns[18] > 1) { $s16ncprobe++;}
    if ($columns[19] > 1) { $s17ncprobe++;}
    if ($columns[20] > 1) { $s18ncprobe++;}
    if ($columns[21] > 1) { $s19ncprobe++;}
    if ($columns[22] > 1) { $s20ncprobe++;}
    if ($columns[23] > 1) { $s21ncprobe++;}
    if ($columns[24] > 1) { $s22ncprobe++;}
    if ($columns[25] > 1) { $s23ncprobe++;}
    if ($columns[26] > 1) { $s24ncprobe++;}
    if ($columns[27] > 1) { $s25ncprobe++;}
    if ($columns[28] > 1) { $s26ncprobe++;}
    if ($columns[29] > 1) { $s27ncprobe++;}
    if ($columns[30] > 1) { $s28ncprobe++;}
    if ($columns[31] > 1) { $s29ncprobe++;}
    if ($columns[32] > 1) { $s30ncprobe++;}
    if ($columns[33] > 1) { $s31ncprobe++;}
    if ($columns[34] > 1) { $s32ncprobe++;}
    if ($columns[35] > 1) { $s33ncprobe++;}
    if ($columns[36] > 1) { $s34ncprobe++;}
    if ($columns[37] > 1) { $s35ncprobe++;}
    if ($columns[38] > 1) { $s36ncprobe++;}
    if ($columns[39] > 1) { $s37ncprobe++;}
    if ($columns[40] > 1) { $s38ncprobe++;}
    if ($columns[41] > 1) { $s39ncprobe++;}
    if ($columns[42] > 1) { $s40ncprobe++;}
        if ($columns[43] > 1) { $s41ncprobe++;}
        if ($columns[44] > 1) { $s42ncprobe++;}
        if ($columns[45] > 1) { $s43ncprobe++;}
        if ($columns[46] > 1) { $s44ncprobe++;}
}
   
   
                 if ($biotype=~/^TUCP/ and $probe =~/^Y$/)
   {
    # split the line into two fields                                                                                                    
    my @columns = split(/\t/, $line);

    if ($columns[3] > 1) {$s1tcpprobe++;}
    if ($columns[4] > 1) { $s2tcpprobe++;}
    if ($columns[5] > 1) { $s3tcpprobe++;}
    if ($columns[6] > 1) { $s4tcpprobe++;}
    if ($columns[7] > 1) { $s5tcpprobe++;}
    if ($columns[8] > 1) { $s6tcpprobe++;}
    if ($columns[9] > 1) { $s7tcpprobe++;}
    if ($columns[10] > 1) { $s8tcpprobe++;}
    if ($columns[11] > 1) { $s9tcpprobe++;}
    if ($columns[12] > 1) { $s10tcpprobe++;}
    if ($columns[13] > 1) { $s11tcpprobe++;}
    if ($columns[14] > 1) { $s12tcpprobe++;}
    if ($columns[15] > 1) { $s13tcpprobe++;}
    if ($columns[16] > 1) { $s14tcpprobe++;}
    if ($columns[17] > 1) { $s15tcpprobe++;}
    if ($columns[18] > 1) { $s16tcpprobe++;}
    if ($columns[19] > 1) { $s17tcpprobe++;}
    if ($columns[20] > 1) { $s18tcpprobe++;}
    if ($columns[21] > 1) { $s19tcpprobe++;}
    if ($columns[22] > 1) { $s20tcpprobe++;}
    if ($columns[23] > 1) { $s21tcpprobe++;}
    if ($columns[24] > 1) { $s22tcpprobe++;}
    if ($columns[25] > 1) { $s23tcpprobe++;}
    if ($columns[26] > 1) { $s24tcpprobe++;}
    if ($columns[27] > 1) { $s25tcpprobe++;}
    if ($columns[28] > 1) { $s26tcpprobe++;}
    if ($columns[29] > 1) { $s27tcpprobe++;}
    if ($columns[30] > 1) { $s28tcpprobe++;}
    if ($columns[31] > 1) { $s29tcpprobe++;}
    if ($columns[32] > 1) { $s30tcpprobe++;}
    if ($columns[33] > 1) { $s31tcpprobe++;}
    if ($columns[34] > 1) { $s32tcpprobe++;}
    if ($columns[35] > 1) { $s33tcpprobe++;}
    if ($columns[36] > 1) { $s34tcpprobe++;}
    if ($columns[37] > 1) { $s35tcpprobe++;}
    if ($columns[38] > 1) { $s36tcpprobe++;}
    if ($columns[39] > 1) { $s37tcpprobe++;}
    if ($columns[40] > 1) { $s38tcpprobe++;}
    if ($columns[41] > 1) { $s39tcpprobe++;}
    if ($columns[42] > 1) { $s40tcpprobe++;}
        if ($columns[43] > 1) { $s41tcpprobe++;}
        if ($columns[44] > 1) { $s42tcpprobe++;}
        if ($columns[45] > 1) { $s43tcpprobe++;}
        if ($columns[46] > 1) { $s44tcpprobe++;}
}
   
   
   
}
}
print "TOTAL BREAK-DOWNS\n";
print"\nprotein_coding\t$s1p\t$s2p\t$s3p\t$s4p\t$s5p\t$s6p\t$s7p\t$s8p\t$s9p\t$s10p\t$s11p\t$s12p\t$s13p\t$s14p\t$s15p\t$s16p\t$s17p\t$s18p\t$s19p\t$s20p\t$s21p\t$s22p\t$s23p\t$s24p\t$s25p\t$s26p\t$s27p\t$s28p\t$s29p\t$s30p\t$s31p\t$s32p\t$s33p\t$s34p\t$s35p\t$s36p\t$s37p\t$s38p\t$s39p\t$s40p\t$s41p\t$s42p\t$s43p\t$s44p\n";
print"lncRNA\t$s1\t$s2\t$s3\t$s4\t$s5\t$s6\t$s7\t$s8\t$s9\t$s10\t$s11\t$s12\t$s13\t$s14\t$s15\t$s16\t$s17\t$s18\t$s19\t$s20\t$s21\t$s22\t$s23\t$s24\t$s25\t$s26\t$s27\t$s28\t$s29\t$s30\t$s31\t$s32\t$s33\t$s34\t$s35\t$s36\t$s37\t$s38\t$s39\t$s40\t$s41\t$s42\t$s43\t$s44\n";
print"TUCP\t$s1tcp\t$s2tcp\t$s3tcp\t$s4tcp\t$s5tcp\t$s6tcp\t$s7tcp\t$s8tcp\t$s9tcp\t$s10tcp\t$s11tcp\t$s12tcp\t$s13tcp\t$s14tcp\t$s15tcp\t$s16tcp\t$s17tcp\t$s18tcp\t$s19tcp\t$s20tcp\t$s21tcp\t$s22tcp\t$s23tcp\t$s24tcp\t$s25tcp\t$s26tcp\t$s27tcp\t$s28tcp\t$s29tcp\t$s30tcp\t$s31tcp\t$s32tcp\t$s33tcp\t$s34tcp\t$s35tcp\t$s36tcp\t$s37tcp\t$s38tcp\t$s39tcp\t$s40tcp\t$s41tcp\t$s42tcp\t$s43tcp\t$s44tcp\n";
print"ncRNA\t$s1nc\t$s2nc\t$s3nc\t$s4nc\t$s5nc\t$s6nc\t$s7nc\t$s8nc\t$s9nc\t$s10nc\t$s11nc\t$s12nc\t$s13nc\t$s14nc\t$s15nc\t$s16nc\t$s17nc\t$s18nc\t$s19nc\t$s20nc\t$s21nc\t$s22nc\t$s23nc\t$s24nc\t$s25nc\t$s26nc\t$s27nc\t$s28nc\t$s29nc\t$s30nc\t$s31nc\t$s32nc\t$s33nc\t$s34nc\t$s35nc\t$s36nc\t$s37nc\t$s38nc\t$s39nc\t$s40nc\t$s41nc\t$s42nc\t$s43nc\t$s44nc\n";
print"pseudogene\t$s1ps\t$s2ps\t$s3ps\t$s4ps\t$s5ps\t$s6ps\t$s7ps\t$s8ps\t$s9ps\t$s10ps\t$s11ps\t$s12ps\t$s13ps\t$s14ps\t$s15ps\t$s16ps\t$s17ps\t$s18ps\t$s19ps\t$s20ps\t$s21ps\t$s22ps\t$s23ps\t$s24ps\t$s25ps\t$s26ps\t$s27ps\t$s28ps\t$s29ps\t$s30ps\t$s31ps\t$s32ps\t$s33ps\t$s34ps\t$s35ps\t$s36ps\t$s37ps\t$s38ps\t$s39ps\t$s40ps\t$s41ps\t$s42ps\t$s43ps\t$s44ps\n\n";


print "ON-TARGET BREAK-DOWNS\n";
print"\nprotein_coding\t$s1pprobe\t$s2pprobe\t$s3pprobe\t$s4pprobe\t$s5pprobe\t$s6pprobe\t$s7pprobe\t$s8pprobe\t$s9pprobe\t$s10pprobe\t$s11pprobe\t$s12pprobe\t$s13pprobe\t$s14pprobe\t$s15pprobe\t$s16pprobe\t$s17pprobe\t$s18pprobe\t$s19pprobe\t$s20pprobe\t$s21pprobe\t$s22pprobe\t$s23pprobe\t$s24pprobe\t$s25pprobe\t$s26pprobe\t$s27pprobe\t$s28pprobe\t$s29pprobe\t$s30pprobe\t$s31pprobe\t$s32pprobe\t$s33pprobe\t$s34pprobe\t$s35pprobe\t$s36pprobe\t$s37pprobe\t$s38pprobe\t$s39pprobe\t$s40pprobe\t$s41pprobe\t$s42pprobe\t$s43pprobe\t$s44pprobe\n";
print"lncRNA\t$s1probe\t$s2probe\t$s3probe\t$s4probe\t$s5probe\t$s6probe\t$s7probe\t$s8probe\t$s9probe\t$s10probe\t$s11probe\t$s12probe\t$s13probe\t$s14probe\t$s15probe\t$s16probe\t$s17probe\t$s18probe\t$s19probe\t$s20probe\t$s21probe\t$s22probe\t$s23probe\t$s24probe\t$s25probe\t$s26probe\t$s27probe\t$s28probe\t$s29probe\t$s30probe\t$s31probe\t$s32probe\t$s33probe\t$s34probe\t$s35probe\t$s36probe\t$s37probe\t$s38probe\t$s39probe\t$s40probe\t$s41probe\t$s42probe\t$s43probe\t$s44probe\n";
print"TUCP\t$s1tcpprobe\t$s2tcpprobe\t$s3tcpprobe\t$s4tcpprobe\t$s5tcpprobe\t$s6tcpprobe\t$s7tcpprobe\t$s8tcpprobe\t$s9tcpprobe\t$s10tcpprobe\t$s11tcpprobe\t$s12tcpprobe\t$s13tcpprobe\t$s14tcpprobe\t$s15tcpprobe\t$s16tcpprobe\t$s17tcpprobe\t$s18tcpprobe\t$s19tcpprobe\t$s20tcpprobe\t$s21tcpprobe\t$s22tcpprobe\t$s23tcpprobe\t$s24tcpprobe\t$s25tcpprobe\t$s26tcpprobe\t$s27tcpprobe\t$s28tcpprobe\t$s29tcpprobe\t$s30tcpprobe\t$s31tcpprobe\t$s32tcpprobe\t$s33tcpprobe\t$s34tcpprobe\t$s35tcpprobe\t$s36tcpprobe\t$s37tcpprobe\t$s38tcpprobe\t$s39tcpprobe\t$s40tcpprobe\t$s41tcpprobe\t$s42tcpprobe\t$s43tcpprobe\t$s44tcpprobe\n";
print"ncRNA\t$s1ncprobe\t$s2ncprobe\t$s3ncprobe\t$s4ncprobe\t$s5ncprobe\t$s6ncprobe\t$s7ncprobe\t$s8ncprobe\t$s9ncprobe\t$s10ncprobe\t$s11ncprobe\t$s12ncprobe\t$s13ncprobe\t$s14ncprobe\t$s15ncprobe\t$s16ncprobe\t$s17ncprobe\t$s18ncprobe\t$s19ncprobe\t$s20ncprobe\t$s21ncprobe\t$s22ncprobe\t$s23ncprobe\t$s24ncprobe\t$s25ncprobe\t$s26ncprobe\t$s27ncprobe\t$s28ncprobe\t$s29ncprobe\t$s30ncprobe\t$s31ncprobe\t$s32ncprobe\t$s33ncprobe\t$s34ncprobe\t$s35ncprobe\t$s36ncprobe\t$s37ncprobe\t$s38ncprobe\t$s39ncprobe\t$s40ncprobe\t$s41ncprobe\t$s42ncprobe\t$s43ncprobe\t$s44ncprobe\n";
print"pseudogene\t$s1psprobe\t$s2psprobe\t$s3psprobe\t$s4psprobe\t$s5psprobe\t$s6psprobe\t$s7psprobe\t$s8psprobe\t$s9psprobe\t$s10psprobe\t$s11psprobe\t$s12psprobe\t$s13psprobe\t$s14psprobe\t$s15psprobe\t$s16psprobe\t$s17psprobe\t$s18psprobe\t$s19psprobe\t$s20psprobe\t$s21psprobe\t$s22psprobe\t$s23psprobe\t$s24psprobe\t$s25psprobe\t$s26psprobe\t$s27psprobe\t$s28psprobe\t$s29psprobe\t$s30psprobe\t$s31psprobe\t$s32psprobe\t$s33psprobe\t$s34psprobe\t$s35psprobe\t$s36psprobe\t$s37psprobe\t$s38psprobe\t$s39psprobe\t$s40psprobe\t$s41psprobe\t$s42psprobe\t$s43psprobe\t$s44psprobe\n\n";

