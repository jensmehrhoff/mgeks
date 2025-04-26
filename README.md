# M-GEKS approach

A Rapid and Efficient Compilation Procedure for the GEKS Index: The Matrix-Based (M-GEKS) Approach

## Purpose of this repo

This repo stores the code to execute the M-GEKS method and also an example benchmark of this method against other datasets. Presented at the [2025 Meeting of the Group of Experts on Consumer Price Indices](https://unece.org/info/events/event/394156)

## Data to support benchmark

The package utilizes analysis-ready data from Dominick's Finer Foods data set. See [Jens Mehrhoff (2018) Promoting the Use of a Publicly Available Scanner Data Set](https://github.com/eurostat/dff/blob/master/docs/dff.pdf) and [dff repository for more info on the data and example code](https://github.com/eurostat/dff)

## Repo structure

```         
├── data
│   ├── raw            <- Folder to place original raw Dominicks's data
│   └── processed      <- Folder to store analysis ready data
│
├── docs               <- Overview of the package and method
│
└── src                <- Source code to compute M-GEKS, run benchmark tests, and code to transform raw to processed data
```