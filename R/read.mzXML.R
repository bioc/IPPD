new.mzXML <- function () 
{
    object = list(header = NULL, parentFile = NULL, dataProcessing = NULL, 
      msInstrument = NULL, separation = NULL, spotting = NULL, 
        scan = vector(mode = "list"))
    class(object) <- "mzXML"
    return(object)
}

read.mzXML <- function(filename) 
{
    Paste = function(...) paste(..., sep = "", collapse = "")
    strtrunc = function(Str, Sub) {
        lp = attr(regexpr(paste(".*", Sub, sep = ""), Str), "match.length")
        return(substring(Str, 1, lp))
      }
    fregexpr = function(pattern, filename) {
        buf.size = 1024
        n = file.info(filename)$size
        pos = NULL
        fp = file(filename, "rb")
        for (d in seq(1, n, by = buf.size)) {
            m = if (n - d > buf.size) 
                buf.size
            else n - d
            p = gregexpr(pattern, readChar(fp, m))[[1]]
            if (p[1] > 0) 
                pos = c(pos, p + d - 1)
        }
        close(fp)
        if (is.null(pos)) 
            pos = -1
        return(pos)
    }
    mzXMLhandlers <- function() {
        obj = new.mzXML()
        iScan = 0
        ParentID = vector(mode = "integer")
        sha1 = vector(mode = "list", length = 2)
        sha1[1] <- sha1[2] <- 0
        OptScanAttr = c("polarity", "scanType", "centroided", 
            "deisotoped", "chargeDeconvoluted", "retentionTime", 
            "ionisationEnergy", "collisionEnergy", "cidGasPressure", 
            "totIonCurrent")
        ToString = function(x, indent = 1) {
            if (is.null(x)) 
                return(NULL)
            spaces = if (indent > 0) 
                Paste(rep("  ", indent))
            else ""
            Name = xmlName(x, TRUE)
            val = xmlValue(x)
            if (Name == "text") 
                return(Paste(spaces, val, "\n"))
            if (!is.null(xmlAttrs(x))) {
                att = paste(names(xmlAttrs(x)), paste("\"", xmlAttrs(x), 
                  "\"", sep = ""), sep = "=", collapse = " ")
                att = paste(" ", att, sep = "")
            }
            else att = ""
            chl = ""
            for (i in xmlChildren(x)) chl = Paste(chl, ToString(i, 
                indent + 1))
            if (chl == "") 
                Str = Paste(spaces, "<", Name, att, "/>\n")
            else Str = Paste(spaces, "<", Name, att, ">\n", chl, 
                spaces, "</", Name, ">\n")
            return(Str)
        }
        CatNodes = function(x, Name, indent = 2) {
            Str = NULL
            for (y in xmlElementsByTagName(x, Name)) Str = paste(Str, 
                ToString(y, indent), sep = "")
            return(Str)
        }
        read.mzXML.scan = function(x) {
            if (is.null(x)) 
                return(NULL)
            if (xmlName(x) != "scan") 
                return(NULL)
            scanOrigin <- precursorMz <- nameValue <- maldi <- mass <- peaks <- NULL
            att = xmlAttrs(x)
            num = as.integer(att["num"])
            msLevel = as.integer(att["msLevel"])
            peaksCount = as.integer(att["peaksCount"])
            msk = names(att) %in% OptScanAttr
            if (sum(msk) == 0) 
                scanAttr = ""
            else {
                scanAttr = paste(names(att[msk]), paste("\"", 
                  att[msk], "\"", sep = ""), sep = "=", collapse = " ")
                scanAttr = paste(" ", scanAttr, sep = "")
            }
            maldi = ToString(x[["maldi"]])
            scanOrigin = CatNodes(x, "scanOrigin", 3)
            nameValue = CatNodes(x, "nameValue", 3)
            precursorMz = CatNodes(x, "precursorMz", 3)
            precursorMz = gsub("\n      ", " ", precursorMz)
            for (y in xmlElementsByTagName(x, "scan")) ParentID[as.integer(xmlAttrs(y)["num"])] <<- num
            y = x[["peaks"]]
            att = xmlAttrs(y)
            peaks = xmlValue(y)
            precision = att["precision"]
            byteOrder = att["byteOrder"]
            pairOrder = att["pairOrder"]
            endian = if (byteOrder == "network") 
                "big"
            else "little"
            if (precision == "32") 
                size = 4
            else if (precision == "64") 
                size = 8
            else stop("read.mzXML.scan: incorrect precision attribute of peaks field")
            if (pairOrder != "m/z-int") 
                warning("read.mzXML.scan: incorrect pairOrder attribute of peaks field")
            if (peaksCount > 0) {
                p = base64decode(peaks, "double", endian = endian, 
                  size = size)
                np = length(p)%/%2
                if (np != peaksCount) 
                  warning("read.mzXML.scan: incorrect 'peakCount' attribute of 'peaks' field: expected ", 
                    peaksCount, ", found ", np, "  ", (3 * ((nchar(peaks) * 
                      size)/4))/2, " (number of scans ", num, ")")
                dim(p) = c(2, np)
                mass = p[1, ]
                peaks = p[2, ]
            }
            return(list(mass = mass, peaks = peaks, num = num, 
                parentNum = num, msLevel = msLevel, scanAttr = scanAttr, 
                maldi = maldi, scanOrigin = scanOrigin, precursorMz = precursorMz, 
                nameValue = nameValue))
          }
        list(mzXML = function(x, ...) {
            y = x[["sha1"]]
            sha1[1] <<- if (!is.null(y)) xmlValue(y) else 0
            x$children = NULL
            obj$header <<- toString(x, terminate = FALSE)
            NULL
        }, msRun = function(x, ...) {
            y = x[["sha1"]]
            sha1[2] <<- if (!is.null(y)) xmlValue(y) else 0
            obj$msInstrument <<- ToString(x[["msInstrument"]], 
                2)
            obj$separation <<- ToString(x[["separation"]], 2)
            obj$spotting <<- ToString(x[["spotting"]], 2)
            obj$parentFile <<- CatNodes(x, "parentFile")
            obj$dataProcessing <<- CatNodes(x, "dataProcessing")
            NULL
        }, scan = function(x, ...) {
            iScan <<- iScan + 1
            obj$scan[[iScan]] <<- read.mzXML.scan(x)
            x$children = NULL
            x
        }, data = function() {
            if (is.null(obj$header)) NULL else list(mzXML = obj, 
                ParentID = ParentID, sha1 = sha1)
        })
    }
    if (!is.character(filename)) 
        stop("read.mzXML: 'filename' has to be a string")
    if (length(filename) > 1) 
        filename = paste(filename, collapse = "")
    sha1File = digest(filename, algo = "sha1", file = TRUE)
    x = xmlTreeParse(file = filename, handlers = mzXMLhandlers(), 
        addAttributeNamespaces = TRUE)$data()
    if (is.null(x)) 
        stop("read.mzXML: This is not mzXML file")
    mzXML = x$mzXML
    sha1Read = x$sha1
    n = length(mzXML$scan)
    NumID = integer(n)
    for (i in 1:n) {
        NumID[i] = mzXML$scan[[i]]$num
        mzXML$scan[[i]]$scanOrigin = paste("<scanOrigin parentFileID='", 
            sha1File, "' num='", NumID[i], "'/>\n", sep = "")
    }
    mzXML$scan = mzXML$scan[order(NumID)]
    for (i in 1:n) if (!is.na(x$ParentID[i])) 
        mzXML$scan[[i]]$parentNum = x$ParentID[i]
    else x$ParentID[i] = mzXML$scan[[i]]$parentNum
    mzXML$scan = mzXML$scan[order(x$ParentID)]
    n = sum(as.integer(lapply(sha1Read, is.character)))
    if (n > 0) {
        if (is.null(sha1Read[[1]])) 
            sha1Read[[1]] = sha1Read[[2]]
        sha1Pos = fregexpr("<sha1>", filename) + 6
        for (i in n) {
            sha1Calc = digest(filename, algo = "sha1", file = TRUE, 
                length = sha1Pos[i] - 1)
            if (sha1Read[[i]] != sha1Calc) 
                warning("Stored and calculated Sha-1 sums do not match (stored '", 
                  sha1Read[[i]], "'; calculated '", sha1Calc, 
                  "')")
        }
    }
    mzXML$header = gsub("/>", ">\n", mzXML$header)
    mzXML$header = gsub("^ +", "", mzXML$header)
    mzXML$header = gsub("[\u0093\u0094]", "\"", mzXML$header)
    mzXML$parentFile = Paste(mzXML$parentFile, "    <parentFile filename='file://", 
        filename, "' fileType='processedData' fileSha1='", sha1File, 
        "'/>\n")
    return(mzXML)
      }
 
