function [H] = fdnEmptyHMatrix(numChr)
    H.s100kb.oe = cell(numChr,1);
    H.s100kb.kr = cell(numChr,1);
    H.s100kb.oeTrim = cell(numChr,1);
    H.s100kb.krTrim = cell(numChr,1);
    H.s100kb.genesTrim = cell(numChr,1);
    H.s100kb.oeTrimBadLocs = cell(numChr,1);
    H.s100kb.centLocStart = cell(numChr,1);
    
    H.s1mb.oe = cell(numChr,1);
    H.s1mb.kr = cell(numChr,1);
    H.s1mb.oeTrim = cell(numChr,1);
    H.s1mb.krTrim = cell(numChr,1);
    H.s1mb.oeTrimBadLocs = cell(numChr,1);
    H.s1mb.chrStart = cell(numChr,1);
    H.s1mb.chrStartTrim = cell(numChr,1);
end

