function format_seq(seq) {
    var lines = seq.split("\n");
    let res = "";
    for (let i = 0; i < lines.length; i+=2) {
        res += lines[i] + "\n";
        if (i+1 < lines.length)
            res += (lines[i+1].match(/.{1,50}/g) || [] ).join("\n") + "\n";
    }
    return res;
}

const copyToClipboard = str => {
  const el = document.createElement('textarea');  // Create a <textarea> element
  el.value = str;                                 // Set its value to the string that you want copied
  el.setAttribute('readonly', '');                // Make it readonly to be tamper-proof
  el.style.position = 'absolute';
  el.style.left = '-9999px';                      // Move outside the screen to make it invisible
  document.body.appendChild(el);                  // Append the <textarea> element to the HTML document
  const selected =
    document.getSelection().rangeCount > 0        // Check if there is any content selected previously
      ? document.getSelection().getRangeAt(0)     // Store selection if found
      : false;                                    // Mark as false to know no selection existed before
  el.select();                                    // Select the <textarea> content
  document.execCommand('copy');                   // Copy - only works as a result of a user action (e.g. click events)
  document.body.removeChild(el);                  // Remove the <textarea> element
  if (selected) {                                 // If a selection existed before copying
    document.getSelection().removeAllRanges();    // Unselect everything on the HTML document
    document.getSelection().addRange(selected);   // Restore the original selection
  }
};

var mapo = {
            'A': 'A',
            'T': 'T',
            'G': 'G',
            'C': 'C',
            'W': 'AT',
            'R': 'AG',
            'K': 'TG',
            'D': 'ATG',
            'M': 'AC',
            'Y': 'CT',
            'H': 'ATC',
            'S': 'CG',
            'V': 'AGC',
            'B': 'CGT',
            'N': 'ACGT'};
var complmap = {
            'A': 'T',
            'T': 'A',
            'G': 'C',
            'C': 'G',
            'W': 'TA',
            'R': 'TC',
            'K': 'AC',
            'D': 'TAC',
            'M': 'TG',
            'Y': 'GA',
            'H': 'TAG',
            'S': 'GC',
            'V': 'TCG',
            'B': 'GCA',
            'N': 'ACGT'};



function compare(str, offset, motif, complementary) {
    for (let i = 0; i < motif.length; i++) {
        let str_idx = complementary ? offset + (motif.length-1-i) : i+offset
        let c = str[str_idx];
        if ("ACGT".indexOf(c) == -1) {
            //console.log("bad string ", str.substring(offset, motif.length+offset))
            return false
        }
        let cc = motif[i];
        if (!complementary && c === cc) continue;
        if ("ATGCWRKDMYHSVBN".indexOf(cc) == -1) {
            return false;
        }
        let s = complementary ? complmap[cc] : mapo[cc];
        //console.log(c, cc);
        if (s.indexOf(c) == -1) {
            return false
        }
    }
    //console.log("good string ", str.substring(offset, motif.length+offset), motif)
    return true
}


function myIndexOf(str, motif, offset, complementary) {
    let fwdIdx = -1;
    for (let i = offset; i < str.length - motif.length + 1; i++) {
        if (compare(str, i, motif, false)) {
            //console.log("good string ", str.substring(i, motif.length+i), motif)
            //console.log("found forward");
            fwdIdx = i;
            break;
        }
    }

    let complIdx = -1;
    //console.log("l " + str + " " + (str.length - motif.length + 1) + " " + motif.length);
    if (complementary)
        for (let i = offset; i < str.length - motif.length + 1; i++) {
            if (compare(str, i, motif, true)) {
                //console.log("found compl " + i);
                complIdx = i;
                break;
            }
        }
    if (fwdIdx == -1)
        return complIdx;
    if (complIdx == -1)
        return fwdIdx;

    let idx = fwdIdx < complIdx ? fwdIdx : complIdx;
    return idx;
}

function split_text(text) {
    let sequences = text.split(">").filter(s => s);

    sequences = sequences.map(s => "> " + s.trim())
    sequences = sequences.map(s => {
        var lines = s.split("\n");
        return [lines[0], lines.splice(1).join("")]
    });
    return sequences
}

function calc_frequencies(text) {
    let splitted = split_text(text);
    let counts = {'A':0, 'T':0, 'G':0, 'C':0};

    let total = 0;
    for (let i = 0; i < splitted.length; i++) {
        let name = splitted[i][0];
        let seq = splitted[i][1];
        total += seq.length;
        seq.split('').forEach(function(s) {
            counts[s] ? counts[s]++ : counts[s] = 1;
        });
    }
    return counts;
}

function count_total_len(text) {
    let splitted = split_text(text);
    let total = 0;
    for (let i = 0; i < splitted.length; i++)
        total += splitted[i][1].trim().length;
    return total;
}

function calc_ratios(text) {
    let counts = calc_frequencies(text);
    let total = count_total_len(text);

    for (var l in counts) {
        counts[l] /= total;
    }
    return counts
}

function prob_in_pos(motif, ratios) {
    var prob1 = 1.0;
    motif.toUpperCase().split("").forEach(function(let) {
        let symbols = mapo[let];
        let curP = 0.0;
        symbols.split("").forEach(function(symb) {
            curP += ratios[symb];
        });
        prob1 *= curP;
    });
    return prob1;
}

function prob_in_sec(prob1, seq_len) {
    var prob2 =  1.0 - Math.exp(-prob1 * (seq_len-7)*2);
    return prob2;
}

function calc_chi2(prob_pos, matches, seq_count) {
    var expected = Math.floor(prob_pos * seq_count);
    if (expected == 0)
        expected = 1
    var x = expected - matches;
    var chi2 = Math.floor(x * x / (expected));
    return chi2;
}

function calc_chi2_double(prob_pos, matches, seq_count) {
    var expected = prob_pos * seq_count;
    var x = expected - matches;
    var chi2 = x * x / (expected);
    return chi2;
}

function unsplit_text(text) {
    let res = "";
    for (let i = 0; i < text.length; i++) {
        let s = text[i][0] + "\n" + text[i][1] + "\n";
        res += s;
    }
    return res;
}

function reverseString(str) {
    // Step 1. Use the split() method to return a new array
    var splitString = str.split(""); // var splitString = "hello".split("");
    // ["h", "e", "l", "l", "o"]

    // Step 2. Use the reverse() method to reverse the new created array
    var reverseArray = splitString.reverse(); // var reverseArray = ["h", "e", "l", "l", "o"].reverse();
    // ["o", "l", "l", "e", "h"]

    // Step 3. Use the join() method to join all elements of the array into a string
    var joinArray = reverseArray.join(""); // var joinArray = ["o", "l", "l", "e", "h"].join("");
    // "olleh"

    //Step 4. Return the reversed string
    return joinArray; // "olleh"
}

function motif_var(motif, complementary) {
    let res = "";

    motif.toUpperCase().split('').forEach(function(s) {
        let possible = (complementary ? complmap[s] : mapo[s]);
        res += possible.charAt(Math.floor(Math.random() * possible.length));
    });


    return complementary? reverseString(res) : res;
}

function inject_motif(text, motif, occurence, complementary) {
    let splitted = split_text(text);

    for (let i = 0; i < splitted.length * occurence; i++) {
        let use_compl = false;
        if (complementary)
            use_compl = Math.floor(Math.random() * 2);
        let mot = motif_var(motif, use_compl);

        let start_pos = Math.floor(Math.random() * (splitted[i][1].length - mot.length+1));
        let sub = splitted[i][1].substring(start_pos, mot.length+start_pos);
        splitted[i][1] = splitted[i][1].replace(sub, mot);
    }
    return unsplit_text(splitted);
}

function gen_seq(num, count, motif, occurence, complementary) {
    var text = ""
    var possible = "ACGT"
    for (var seq = 0; seq < num; seq++) {
        text += "> " + seq + "\n";
        for (var i = 0; i < count; i++) {
            text += possible.charAt(Math.floor(Math.random() * possible.length));
        }
        text += "\n"
    }
    text += "\n"

    if (motif) {
        text = inject_motif(text, motif, occurence, complementary);
    }

    return format_seq(text);
}
