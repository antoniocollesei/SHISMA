CUSTOM_CONFIG_A3 = [
    {
        "window_size": 2,
        "stride": 1,
        "dilation": 1,
        "word_length": 2,
        "alphabet_size": 3,
    },
    {
        "window_size": 2,
        "stride": 1,
        "dilation": 2,
        "word_length": 2,
        "alphabet_size": 3,
    },
    {
        "window_size": 3,
        "stride": 1,
        "dilation": 1,
        "word_length": 3,
        "alphabet_size": 3,
    },
    {
        "window_size": 3,
        "stride": 1,
        "dilation": 2,
        "word_length": 3,
        "alphabet_size": 3,
    },
    {
        "window_size": 4,
        "stride": 1,
        "dilation": 1,
        "word_length": 2,
        "alphabet_size": 3,
    },
    {
        "window_size": 4,
        "stride": 1,
        "dilation": 1,
        "word_length": 4,
        "alphabet_size": 3,
    },
    {
        "window_size": 5,
        "stride": 1,
        "dilation": 1,
        "word_length": 5,
        "alphabet_size": 3,
    },
    {
        "window_size": 6,
        "stride": 1,
        "dilation": 1,
        "word_length": 2,
        "alphabet_size": 3,
    },
    {
        "window_size": 6,
        "stride": 1,
        "dilation": 1,
        "word_length": 3,
        "alphabet_size": 3,
    },
    {
        "window_size": 6,
        "stride": 1,
        "dilation": 1,
        "word_length": 6,
        "alphabet_size": 3,
    },
]

CUSTOM_CONFIG_A3_NO_DILATION = [i for i in CUSTOM_CONFIG_A3 if i["dilation"] == 1]

CUSTOM_CONFIG_A3_NO_DILATION_WINDOW_SIZE_2_3_4 = [
    i for i in CUSTOM_CONFIG_A3_NO_DILATION if i["window_size"] in [2, 3, 4]
]
