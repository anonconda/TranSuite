import sys


def get_max_atg(orf_data, locus_model):

    max_orf_length = max(orf_data[1])

    max_orf_models, max_orf_atg = ([] for _ in range(2))
    for i, orf_length in enumerate(orf_data[1]):

        if orf_length == max_orf_length:
            transcript_id = orf_data[0][i]
            max_orf_models.append(transcript_id)

            atg = locus_model.transcript_dict[transcript_id].get_atg()
            max_orf_atg.append(atg)

    if not len(set(max_orf_atg)) == 1:
        # If two ATG gave the same longest ORF, choose the closest to the start of the transcript
        transcript_ID, transcript_model = sorted(locus_model.transcript_dict.items())[0]

        if transcript_model.sense == "+":
            max_orf_atg = [min(max_orf_atg)]
        elif transcript_model.sense == "-":
            max_orf_atg = [max(max_orf_atg)]
        else:
            sys.exit(f"Transcript {transcript_ID} strand must be either + or -, not {transcript_model.sense}")

    return max_orf_atg[0], max_orf_models


def check_atg_universal(atg, locus_model):

    transcript_dict = locus_model.transcript_dict
    for transcript_id in transcript_dict:
        atg_present = False

        transcript_model = transcript_dict[transcript_id]
        exon_list = transcript_model.exon_list

        for exon in exon_list:
            if atg in range(exon[0], exon[1] + 1):
                atg_present = True

        if atg_present is True:
            return True

    # Currently this method is disabled (by always returning True)
    return True


def update_orf_data(orf_data):

    # Remove all instances of current ATG
    new_orf_data = [[], []]

    try:
        max_orf_length = max(orf_data[1])
    except:
        return 0

    for i, orf_length in enumerate(orf_data[1]):
        if orf_length != max_orf_length:
            new_orf_data[0].append(orf_data[0][i])
            new_orf_data[1].append(orf_data[1][i])

    orf_data = new_orf_data

    return orf_data


def characterise_max_orfs(locus_dict):

    ambiguous_locus, max_orf_transcripts, max_orf_counts = ([] for _ in range(3))
    for locus_id in locus_dict:
        locus_model = locus_dict[locus_id]
        orf_data = locus_model.get_orf_lengths()

        atg, max_orf_models = get_max_atg(orf_data, locus_model)

        locus_model.rep_atg = atg
        locus_dict.update({locus_id: locus_model})

        max_orf_transcripts.extend(max_orf_models)
        max_orf_counts.append(len(max_orf_models))

    return max_orf_transcripts, max_orf_counts, locus_dict


def remove_max_orf(locus_dict):

    for locus_id in list(locus_dict.keys()):
        locus_model = locus_dict[locus_id]

        orf_data = locus_model.get_orf_lengths()
        orf_data = update_orf_data(orf_data)

        locus_model.transcript_dict.orf

        try:
            locus_dict.update({locus_id: locus_model})
        except:
            del locus_dict[locus_id]

    return locus_dict
