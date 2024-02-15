function E = SourceFct(t,InputParas)

if isfield(InputParas,'rep')        %If structure InputParas has a field named rep
    n = floor(t/InputParas.rep);
    t = t-n*InputParas.rep;
end

if ~isstruct(InputParas)            %Passes if not a structure
    E = InputParas;
else                                %A structure without rep as a field
    E = InputParas.E0*exp(-(t-InputParas.t0)^2/InputParas.wg^2)*exp(1i*(InputParas.we*t + InputParas.phi));
end

end
