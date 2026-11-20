function k = prt_require_transverse(E,k)
% Check global complex field against real ray direction, then normalize k.
assert(size(E,1)==3 && isequal(size(E),size(k)) && isreal(k) && ...
    all(isfinite(k),'all') && all(isfinite(E),'all'),'三维电场和方向须为有限数值。');
lengths=vecnorm(k);
assert(all(lengths>0 & isfinite(lengths)),'光线方向不能为零。');
k=k./lengths;
bad=find(abs(sum(k.*E,1))>1e-9*vecnorm(E),1);
assert(isempty(bad),'第 %d 个三维采样的电场不垂直于光线方向：须满足 k·E=0。',bad);
end
