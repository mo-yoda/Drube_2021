import pandas as pd


def main():
    path = "C:/path/to/folder/"
    # path_laptop = "C:/path/to/folder/"

    export_csv(path_to_folder=path)
    # export_csv(path_to_folder=path_laptop)

def ident_p(path_to_folder):
    # input: table with GPCR, barr, GRK group, barr class, IL3 or C term, AA sequence
    data = pd.read_csv(path_to_folder + "C_termini_AND_IL3.csv", sep=";")
    print(data)

    # add column for length of AA sequence
    length = []
    for row in data["seq"]:
        length.append(len(row))
    data["length"] = length

    # count P (Ser/Thr) and N (Asp/Glu) sites
    p_sites = []
    n_sites = []
    for seq in data["seq"]:
        p_sites.append(seq.count("S") + seq.count("T"))
        n_sites.append(seq.count("E") + seq.count("D"))
    data["p_count"] = p_sites
    data["n_count"] = n_sites

    # index - row of sequence, to identify GPCR and structure (IL3 or Cterm)
    # pos - (position) position in the respective structure

    # single P positions
    only_p = ["S", "T"]
    pos_P = []
    indices_P = []

    for row, seq in enumerate(data["seq"]):
        for pos, AA in enumerate(seq):
            if AA in only_p:
                pos_P.append(pos)
                indices_P.append(row)

    p_site = ["S", "T", "E", "D"]
    print(data.columns)

    indices_PPP = []
    pos_PPP = []
    indices_PXPP = []
    pos_PXPP = []
    indices_PXPXXP = []
    pos_PXPXXP = []
    indices_PXXPXXP = []
    pos_PXXPXXP = []

    # pattern recognition
    for row, seq in enumerate(data["seq"]):

        for pos, AA in enumerate(seq):
            if AA in p_site:
                # P start
                try:
                    if seq[pos + 1] in p_site:
                        # print("PP")
                        try:
                            if seq[pos + 2] in p_site:
                                # print("PPP")
                                indices_PPP.append(row)
                                pos_PPP.append(pos)
                            elif seq[pos + 3] in p_site:
                                # print("PPXP")
                                indices_PXPP.append(row)
                                pos_PXPP.append(pos)
                            elif seq[pos + 2] in p_site:
                                # print("PXP")
                                if seq[pos + 3] in p_site:
                                    # print("PXPP")
                                    indices_PXPP.append(row)
                                    pos_PXPP.append(pos)
                                elif seq[pos + 5] in p_site:
                                    # print("PXPXXP")
                                    indices_PXPXXP.append(row)
                                    pos_PXPXXP.append(pos)
                        except IndexError:
                            continue
                    else:
                        try:
                            if seq[pos + 3] in p_site:
                                # print("PXXP")
                                if seq[pos + 5] in p_site:
                                    # print("PXXPXP")
                                    indices_PXPXXP.append(row)
                                    pos_PXPXXP.append(pos)
                                if seq[pos + 6] in p_site:
                                    # print("PXXPXXP")
                                    indices_PXXPXXP.append(row)
                                    pos_PXXPXXP.append(pos)
                        except IndexError:
                            continue
                except IndexError:
                    continue

    indices_five_in_eight = []
    pos_five_in_eight = []
    # moving frame 5x P in 8 AA
    for row, seq in enumerate(data["seq"]):
        for pos, AA in enumerate(seq):
            try:
                seq[pos + 8]
            except IndexError:
                pass
            temp_seq = seq[pos: pos + 8]
            p_count = temp_seq.count("S") + temp_seq.count("T") + temp_seq.count("N")
            if p_count >= 5:
                # more than 5x P in 8 AA
                indices_five_in_eight.append(row)
                pos_five_in_eight.append(pos)


    # position table for single P (Ser or Thr)
    GPCR = []
    barr = []
    GRK = []
    arr_class = []
    IL3_or_Cterm = []
    length = []
    for index in indices_P:
        GPCR.append(data["GPCR"][index])
        barr.append(data["barr"][index])
        GRK.append(data["GRK"][index])
        arr_class.append(data["arr_class"][index])
        IL3_or_Cterm.append(data["IL3_or_Cterm"][index])
        length.append(data["length"][index])
    P_positions = pd.DataFrame()
    P_positions["GPCR"] = GPCR
    P_positions["barr"] = barr
    P_positions["GRK"] = GRK
    P_positions["arr_class"] = arr_class
    P_positions["IL3_or_Cterm"] = IL3_or_Cterm
    P_positions["length"] = length
    P_positions["PPP_position"] = pos_P

    # position table for PPP
    GPCR = []
    barr = []
    GRK = []
    arr_class = []
    IL3_or_Cterm = []
    length = []
    for index in indices_PPP:
        GPCR.append(data["GPCR"][index])
        barr.append(data["barr"][index])
        GRK.append(data["GRK"][index])
        arr_class.append(data["arr_class"][index])
        IL3_or_Cterm.append(data["IL3_or_Cterm"][index])
        length.append(data["length"][index])
    PPP_positions = pd.DataFrame()
    PPP_positions["GPCR"] = GPCR
    PPP_positions["barr"] = barr
    PPP_positions["GRK"] = GRK
    PPP_positions["arr_class"] = arr_class
    PPP_positions["IL3_or_Cterm"] = IL3_or_Cterm
    PPP_positions["length"] = length
    PPP_positions["PPP_position"] = pos_PPP

    # position table for PXPP
    GPCR = []
    barr = []
    GRK = []
    arr_class = []
    IL3_or_Cterm = []
    length = []
    for index in indices_PXPP:
        GPCR.append(data["GPCR"][index])
        barr.append(data["barr"][index])
        GRK.append(data["GRK"][index])
        arr_class.append(data["arr_class"][index])
        IL3_or_Cterm.append(data["IL3_or_Cterm"][index])
        length.append(data["length"][index])
    PXPP_positions = pd.DataFrame()
    PXPP_positions["GPCR"] = GPCR
    PXPP_positions["barr"] = barr
    PXPP_positions["GRK"] = GRK
    PXPP_positions["arr_class"] = arr_class
    PXPP_positions["IL3_or_Cterm"] = IL3_or_Cterm
    PXPP_positions["length"] = length
    PXPP_positions["PXPP_position"] = pos_PXPP

    # position table for PXPXXP
    GPCR = []
    barr = []
    GRK = []
    arr_class = []
    IL3_or_Cterm = []
    length = []
    for index in indices_PXPXXP:
        GPCR.append(data["GPCR"][index])
        barr.append(data["barr"][index])
        GRK.append(data["GRK"][index])
        arr_class.append(data["arr_class"][index])
        IL3_or_Cterm.append(data["IL3_or_Cterm"][index])
        length.append(data["length"][index])

    PXPXXP_positions = pd.DataFrame()
    PXPXXP_positions["GPCR"] = GPCR
    PXPXXP_positions["barr"] = barr
    PXPXXP_positions["GRK"] = GRK
    PXPXXP_positions["arr_class"] = arr_class
    PXPXXP_positions["IL3_or_Cterm"] = IL3_or_Cterm
    PXPXXP_positions["length"] = length
    PXPXXP_positions["PXPXXP_position"] = pos_PXPXXP

    # position table for PXXPXXP
    GPCR = []
    barr = []
    GRK = []
    arr_class = []
    IL3_or_Cterm = []
    length = []
    for index in indices_PXXPXXP:
        GPCR.append(data["GPCR"][index])
        barr.append(data["barr"][index])
        GRK.append(data["GRK"][index])
        arr_class.append(data["arr_class"][index])
        IL3_or_Cterm.append(data["IL3_or_Cterm"][index])
        length.append(data["length"][index])

    PXXPXXP_positions = pd.DataFrame()
    PXXPXXP_positions["GPCR"] = GPCR
    PXXPXXP_positions["barr"] = barr
    PXXPXXP_positions["GRK"] = GRK
    PXXPXXP_positions["arr_class"] = arr_class
    PXXPXXP_positions["IL3_or_Cterm"] = IL3_or_Cterm
    PXXPXXP_positions["length"] = length
    PXXPXXP_positions["PXXPXXP_position"] = pos_PXXPXXP

    # position table for five P in 8 AA
    GPCR = []
    barr = []
    GRK = []
    arr_class = []
    IL3_or_Cterm = []
    length = []
    for index in indices_five_in_eight:
        GPCR.append(data["GPCR"][index])
        barr.append(data["barr"][index])
        GRK.append(data["GRK"][index])
        arr_class.append(data["arr_class"][index])
        IL3_or_Cterm.append(data["IL3_or_Cterm"][index])
        length.append(data["length"][index])

    fine_positions = pd.DataFrame()
    fine_positions["GPCR"] = GPCR
    fine_positions["barr"] = barr
    fine_positions["GRK"] = GRK
    fine_positions["arr_class"] = arr_class
    fine_positions["IL3_or_Cterm"] = IL3_or_Cterm
    fine_positions["length"] = length
    fine_positions["five_in_eight_position"] = pos_five_in_eight

    # count tables (how many occurrances of P/cluster/pattern)
    PPP_count = []
    PXPP_count = []
    PXPXXP_count = []
    PXXPXXP_count = []
    five_in_eight_count = []
    for index in range(len(data["GPCR"])):
        PPP_count.append(indices_PPP.count(index))
        PXPP_count.append(indices_PXPP.count(index))
        PXPXXP_count.append(indices_PXPXXP.count(index))
        PXXPXXP_count.append(indices_PXXPXXP.count(index))
        five_in_eight_count.append(indices_five_in_eight.count(index))

    data["PPP_count"] = PPP_count
    data["PXPP_count"] = PXPP_count
    data["PXPXXP_count"] = PXPXXP_count
    data["PXXPXXP_count"] = PXXPXXP_count
    data["five_in_eight"] = five_in_eight_count

    return data, PPP_positions, PXPP_positions, PXPXXP_positions, PXXPXXP_positions, fine_positions, P_positions


def export_csv(path_to_folder):
    results = ident_p(path_to_folder)
    results[0].to_csv(path_to_folder + "moving_frame_count.csv", index=False)
    results[1].to_csv(path_to_folder + "PPP_pos.csv", index=False)
    results[2].to_csv(path_to_folder + "PXPP_pos.csv", index=False)
    results[3].to_csv(path_to_folder + "PXPXXP_pos.csv", index=False)
    results[4].to_csv(path_to_folder + "PXXPXXP_pos.csv", index=False)
    results[5].to_csv(path_to_folder + "fine_pos.csv", index=False)
    results[6].to_csv(path_to_folder + "P_pos.csv", index=False)

    results[0].to_excel(path_to_folder + "moving_frame_count.xlsx", index=False)
    results[1].to_excel(path_to_folder + "PPP_pos.xlsx", index=False)
    results[2].to_excel(path_to_folder + "PXPP_pos.xlsx", index=False)
    results[3].to_excel(path_to_folder + "PXPXXP_pos.xlsx", index=False)
    results[4].to_excel(path_to_folder + "PXXPXXP_pos.xlsx", index=False)
    results[5].to_excel(path_to_folder + "fine_pos.xlsx", index=False)
    results[6].to_excel(path_to_folder + "P_pos.xlsx", index=False)


if __name__ == '__main__':
    main()
