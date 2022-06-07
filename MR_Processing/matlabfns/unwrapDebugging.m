clear all; close all;

load('tempphase.mat');

uwPh = WeightedUnwrap3(hCI_phx.*mask,ones(size(hCI_phx)),mask);
%uwPh2 = WeightedUnwrap3(uwPh.*mask,ones(size(hCI_phx)),mask);